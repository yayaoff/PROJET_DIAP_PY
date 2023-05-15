#include "eigen.h"
#include <math.h>
#include <stdlib.h>
#include "lu.h"
#include "elasticity.h"
#include "../../gmsh-sdk/include/gmshc.h"
#include "design.h"

typedef struct{
    int argc;
    int k;
    double r1;
    double r2;
    double e;
    double l; 
    double meshSizeFactor;
    char *filename;
    char **argv;
}param_t;

double *compute_freq(int tag_vis, param_t *param);