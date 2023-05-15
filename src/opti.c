#include "../headers/opti.h"

double *compute_freq(int tag_vis, param_t *param){
    // Number of vibration modes to find
    int k = param->k; 
    double r1 = param->r1 ; double r2 = param->r2; double e = param->e; double l = param->l; 
    double meshSizeFactor = param->meshSizeFactor;

    double *frequencies = malloc(sizeof(double)*k);

    // Define physical constants
    double E = 0.7e11;  // Young's modulus for Aluminum
    double nu = 0.3;    // Poisson coefficient
    double rho = 3000;  // Density of Aluminum

    // Initialize Gmsh and create geometry
    int ierr;
    gmshInitialize(0, NULL, 0, 0, &ierr);

    // designTuningFork(6e-3, 11e-3, 38e-3, 82e-3, 0.3, NULL);
    designTuningFork(r1, r2, e, l, meshSizeFactor, NULL);

    // Assemble the 2 matrices of the linear elasticity problem: 
    // M is the mass matrix
    // K is the stiffness matrix
    Matrix *K, *M;
    size_t* boundary_nodes;
    size_t n_boundary_nodes;
    double * coord;
    assemble_system(&K, &M, &coord, &boundary_nodes, &n_boundary_nodes, E, nu, rho);

    // Remove lines from matrix that are boundary
    Matrix *K_new;
    Matrix *M_new;
    remove_bnd_lines(K, M, boundary_nodes, n_boundary_nodes, &K_new, &M_new, NULL);

    // Compute M := K^{-1} M 
    int n = K_new->n;
    lu(K_new);
    double Mj[n]; // column j
    for (int j=0; j<n; j++){ // column by column
    for (int i=0; i<n; i++) Mj[i] = M_new->a[i][j];
    solve(K_new, Mj);
    for (int i=0; i<n; i++){
        M_new->a[i][j] = Mj[i];
    }
    } // -> Now, in M_new we have stored K^-1M

    // Power iteration + deflation to find k largest eigenvalues
    Matrix * A = M_new; // alias for conciseness
    double * v = malloc(A->m * sizeof(double));
    double lambda, freq;
    FILE * file = fopen(param->argv[2], "w"); // open file to write frequencies
    for(int ki = 0; ki < k; ki++) {
    lambda = power_iteration(A, v);
    freq = 1./(2*M_PI*sqrt(lambda));

    fprintf(file, "%.9lf ", freq);

    frequencies[ki] = freq;

    printf("lambda = %.9e, f = %.3lf\n", lambda, freq);
    // Deflate matrix
    for(int i = 0; i < A->m; i++)
        for(int j = 0; j < A->n; j++)
        A->a[i][j] -= lambda * v[i] * v[j];

    // Put appropriate BC and plot
    double * vall = calloc(K->m, sizeof(double));
    int iv = 0, i_bnd = 0; 
    for(int i = 0; i < K->m/2; i++) {
        if(i_bnd < n_boundary_nodes && i == boundary_nodes[i_bnd]) {
        i_bnd++;
        continue;
        }
        vall[2*(i)]   = v[2*iv];
        vall[2*(i)+1] = v[2*iv+1];
        iv++;
    }
    visualize_in_gmsh(vall, K->m/2);
    }

    fclose(file);

    gmshFltkRun(&ierr);

    free_matrix (K);
    free_matrix (M);
    free_matrix (K_new);
    free_matrix (M_new);
    free(boundary_nodes);

    return frequencies;
}