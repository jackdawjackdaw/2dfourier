#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "decompf2d.h"

/**
 * @file example-driver.c
 * @Author ccs
 * @date 14.09.2012
 * @brief example driver routine, reads an input file and computes the decomposition moments Amn
 * 
 * 
 * reads an input file and computes the Amn and reports the various moments defined.
 * the number of points to use in the grid and the maximum moments to compute are taken as arguments
 * run as: 
 * f2d-driver-example ./example/example-event.dat 200 8 8 
 * 
 * 
 * input file is a *2d* grid of energy density representing the initial-state of an event (or really any part of the evolution)
 * 
 * each line:
 * ix iy edensity
 *
 * ix: (int) x coordinate of the cell
 * iy: (int) y coordinate of the cell
 * edensity: (double) energy density in the cell
 *
 * 
 */ 


int main (int argc, char* argv[]){
  FILE *fptr = NULL;
  char buffer[256];
  gsl_matrix *midplane;
  gsl_matrix *AmnRe, *AmnIm;
  int i, j;
  int mmax=2, nmax=4;
  int npts = 200;
  int xtemp, ytemp;
  double etemp = 0.0;

  double cmx = 0.0, cmy = 0.0;
  double dx = 0.0; 
  double mod = 0.0;
  
  if(argc < 5){
    fprintf(stderr, "# requries args:\n# <s:input-file> <i:npts> <i:nmax> <i:mmax> ");
    exit(-1);
  }

  sprintf(buffer, "%s", argv[1]);
  printf("# reading from: %s\n", buffer);


  npts = atoi(argv[2]);
  assert(npts > 0);
  printf("# input grid size: (%d x %d)\n", npts, npts);

  mmax = atoi(argv[3]);
  assert(mmax > 0 && mmax < MAXMDECOMP);
  nmax = atoi(argv[4]);
  assert(nmax > 0 && nmax < MAXNDECOMP);
  printf("# moments computed up to: %d %d\n", mmax, nmax);
  
  dx = 2*(fabs(xmin))/((double)npts-1);
  /* allocate a 2d gsl matrix that will store the event plane energy density 
   * each cell in this grid should correspond to a hydro cell or similar */
  midplane = gsl_matrix_alloc(npts, npts);

  fptr = fopen(buffer, "r");
  if(!fptr){
    fprintf(stderr, "# cannot open %s\n", buffer);
    exit(-1);
  }

  while(fscanf(fptr, "%d %d %lf",
               &xtemp, &ytemp, &etemp) != EOF){
    gsl_matrix_set(midplane, xtemp - 1, ytemp -1, etemp);
  };
  
  fclose(fptr);

  /* allocate matrices to hold the real and im parts of the coeffs */
  AmnRe = gsl_matrix_alloc((2*mmax+1), nmax);
  AmnIm = gsl_matrix_alloc((2*mmax+1), nmax);

  /* get the cm of this event */
  compute_com(midplane, npts, &cmx, &cmy);

  /* print out the cm in the scaled coors and the indicies of the actual grid cell */
  printf("# com: %lf %lf (%d %d)\n", cmx, cmy, (int)((cmx+(fabs(xmin)))/dx), (int)((cmy+(fabs(xmin)))/dx));

  /* do the decomp  */
  compute_amn(mmax, nmax, midplane, npts, AmnRe, AmnIm, cmx, cmy); // first compute with cm at origin of coords

  /* print out the coeffs and compute the L2 norm */
  for(i = 0; i < (2*mmax+1); i++){
    for(j = 0; j < nmax; j++){
      printf("(%d %d)[%lf %lf] ", -mmax+i, j, gsl_matrix_get(AmnRe, i, j), gsl_matrix_get(AmnIm, i,j));
    }
    printf("\n");
  }
  printf("L2: %lf\n", compute_l2(mmax, nmax, AmnRe, AmnIm));
  printf("M1: %lf\n", compute_m1(mmax, nmax, AmnRe, AmnIm));
  printf("H1: %lf\n", compute_h1(mmax, nmax, AmnRe, AmnIm));
  printf("Rsq: %lf\n", compute_rsq(mmax, nmax, AmnRe, AmnIm));
  
  
  gsl_matrix_free(AmnRe);
  gsl_matrix_free(AmnIm);
  gsl_matrix_free(midplane);
  return EXIT_SUCCESS;
}
