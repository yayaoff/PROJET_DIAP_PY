#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../gmsh-sdk/include/gmshc.h"
#include "headers/matrix.h"
#include "headers/elasticity.h"
#include "math.h"
#include "headers/lu.h"
#include "headers/design.h"
#include "headers/eigen.h"
#include "headers/opti.h"

#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })


int main (int argc, char *argv[]) {

  if (argc < 2){
    printf("Usage: \n"
			"./project <k> <out>\n" 
			"---------------------------- \n\n"
			"- k is the number of frequencies to compute. \n "
			"- out is the output file to write the frequencies. \n "
      "\n");
		return -1;
  } 

  param_t *param = malloc(sizeof(param_t));
  param->k = atoi(argv[1]); param->argc = argc; param->argv = argv; param->filename = argv[2];
  //param->r1 = 6e-3; param->r2 = 11e-3; param->e = 38e-3; param->l = 478e-4; param->meshSizeFactor = 0.3;
  param->r1 = 6e-3; param->r2 = 11e-3; param->e = 38e-3; param->l = 82e-3; param->meshSizeFactor = 0.3;
  // designTuningFork(6e-3, 11e-3, 38e-3, 82e-3, 0.3, NULL);
  double *freq_init=compute_freq(1,param); 
  printf("Initial frequence : %f\n",freq_init[0]);


  return 0;
}
