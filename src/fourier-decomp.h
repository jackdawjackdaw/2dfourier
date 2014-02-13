#ifndef _INC_FDECOMP_
#define _INC_EDECOMP_

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sf.h>

/**
 * ccs, 14.09.2012
 * functions to compute the Amn decomposition for a given array of points
 * 08.10.2012
 * now, how can we use this in R?
 */
void compute_amn(int mmax, int nmax, gsl_matrix *array, int npts, gsl_matrix* AmnReal, gsl_matrix* AmnIm, double cmx, double cmy);
void compute_com(gsl_matrix *array, int npts,double* cmx, double* cmy);
void setup_fdecomp_(int *ngrid);
void free_fdecomp_();
void fill_grid_(int *i, int *j, double *val);
void do_fdecomp_(int *nmax_in, int *mmax_in, char* outname[]);

/* sets the scale of the box*/
double xmin = -1;


/** globals for defining the working grid */
gsl_matrix* grid;
int nptsGrid;


/** defines for sanity checking in do_fdecomp */

#define MAXMDECOMP 24
#define MMAXDEFAULT 8
#define MAXNDECOMP 24
#define NMAXDEFAULT 8


#endif