#include "stdio.h"
#include "stdlib.h"
#include "math.h"
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
void do_fdecomp_();

/* sets the scale of the box
 */
double xmin = -1;

gsl_matrix* grid;
int nptsGrid;

// fns to be called from fortran
/**
 * setup the decomp from fortran, this allocs the variable grid
 */
void setup_fdecomp_(int* ngrid)
{
	nptsGrid = *ngrid;
	printf("#(c) npts:%d\n", nptsGrid);
	grid = gsl_matrix_alloc(nptsGrid, nptsGrid);
	gsl_matrix_set_zero(grid); // clear the grid
}

void free_fdecomp_()
{
	gsl_matrix_free(grid);
}

void fill_grid_(int *i, int *j, double *val)
{
	//printf("#(c) (%d %d) = %lf\n", *i-1, *j-1, *val);
	gsl_matrix_set(grid, *i-1, *j-1, *val);
}

// how can we label the output files?
// should send evt number here
void do_fdecomp_()
{
	int mmax=8;
	int nmax=8;
	int i,j;
	double cmx =0.0, cmy=0.0;
	FILE* fptr;
	gsl_matrix * amnRe = gsl_matrix_alloc(2*mmax+1, nmax);
	gsl_matrix * amnIm = gsl_matrix_alloc(2*mmax+1, nmax);

	compute_com(grid, nptsGrid, &cmx, &cmy);
	
	printf("#(c) cm (%lf %lf)\n", cmx, cmy);
	fprintf(stderr, "# started computing fdecomp\n");
	compute_amn(mmax, nmax, grid, nptsGrid, amnRe, amnIm, cmx, cmy);
	
	// need to name this properly...
	fptr = fopen("fdecomp.dat", "a");
	fprintf(stderr, "# started writing fdecomp\n");
	// dump the decomp to file
	fprintf(fptr, "#\n"); // sep events with a #
	for(i = 0; i < (2*mmax+1); i++){
		for(j = 0; j < nmax; j++){
			fprintf(fptr, "%d %d %lf %lf\n", -mmax+i, j, gsl_matrix_get(amnRe, i, j), gsl_matrix_get(amnIm, i,j));
		}
	}
	fclose(fptr);
	
	fprintf(stderr, "# started computing fdecomp\n");
	
	gsl_matrix_free(amnRe);
	gsl_matrix_free(amnIm);

}


// fns to do the decomp
/**
 * compute the fourier decomposition of array to -m:m in angular components and 1:nmax in radial components
 * 
 * array should a pointer to an initialized npts x npts grid containing the data we want to decompose
 * AmnReal and AmnIm should be allocated arrays of size (2*mmax+1) * nmax which will be filled with the coefficients
 * 
 */ 
void compute_amn(int mmax, int nmax, gsl_matrix *array, int npts, gsl_matrix* AmnReal, gsl_matrix* AmnIm, double cmx, double cmy)
{
	int i,j,k,l;
	int nm, nn;
	int mtemp, ntemp;
	double dx;// = 2/((double)npts-1);
	double dxy;// = 4/pow(((double)npts-1),2.0);
	double xv, yv;
	double coeff = 0;
	double phiMod, phiRe, phiIm;
	double ftemp;
	double AmnRealAcc, AmnImAcc;
	// for compensated summation
	double alphaRe, alphaIm;
	double epsRe, epsIm;

	double rzero = fabs(xmin);
	
	dx = 2*rzero/((double)npts-1);
	dxy = 4*pow(rzero,2.0)/pow(((double)npts-1),2.0);
	//printf("# xmin: %lf dx: %lf dxy: %lf\n", xmin, dx, dxy);
	
	gsl_vector *xvec = gsl_vector_alloc(npts);
	gsl_matrix * rMat = gsl_matrix_alloc(npts, npts);
	gsl_matrix * thMat = gsl_matrix_alloc(npts, npts);
	gsl_matrix *lamMat = NULL;

	gsl_matrix_set_zero(rMat);
	gsl_matrix_set_zero(thMat);
	gsl_vector_set_zero(xvec);

	nm = 2*mmax+1;
	nn = nmax;

	lamMat = gsl_matrix_alloc(nm, nn);

	// fill in r and Theta matrices
	for(i = 0; i < npts; i++)
		gsl_vector_set(xvec ,i, xmin + dx*i);
	
	for(i = 0; i < npts; i++){			
		xv = gsl_vector_get(xvec, i);
		for(j = 0; j < npts;j ++){
			yv = gsl_vector_get(xvec, j);
			gsl_matrix_set(rMat, i, j, sqrt((xv-cmx)*(xv-cmx) + (yv-cmy)*(yv-cmy)));
			gsl_matrix_set(thMat, i, j, atan2((yv-cmy), (xv-cmx)));
		}
	}

	// fill in lambda matrix
	for(i=0; i < nm; i++){
		for(j = 0; j < nn; j++){
			ntemp = j + 1;
			mtemp = -1.0 * mmax + i;
			gsl_matrix_set(lamMat, i, j, gsl_sf_bessel_zero_Jnu(fabs(mtemp), ntemp));
		}
	}
	
	for(i = 0; i < nm; i++){
		for(j = 0; j < nn; j++){
			AmnImAcc = 0.0;
			AmnRealAcc = 0.0;
			epsRe = 0.0;
			epsIm = 0.0;
			
			ntemp = j + 1;
			mtemp = -1.0*mmax + i;
			// note that we have to scale the coeff by rzero, then the system is properly scale invariant
			coeff = pow(rzero,2)*sqrt(M_PI)*gsl_sf_bessel_Jn(fabs(mtemp)+1, gsl_matrix_get(lamMat, i, j));
			
			// now loop over the grid, a lot
			for(k = 0; k < npts; k++){
				for(l = 0; l < npts; l++){
					phiMod = gsl_sf_bessel_Jn(mtemp, gsl_matrix_get(lamMat, i, j)*gsl_matrix_get(rMat, k, l)/rzero) / coeff;
					phiRe = phiMod * cos(mtemp*gsl_matrix_get(thMat, k, l));
					phiIm = phiMod * sin(mtemp*gsl_matrix_get(thMat, k, l));
					ftemp = gsl_matrix_get(array, k, l);

					// compensated summation
					alphaRe = AmnRealAcc;
					epsRe += ftemp * phiRe;
					AmnRealAcc = alphaRe + epsRe;
					epsRe += (alphaRe - AmnRealAcc);

					alphaIm = AmnImAcc;
					epsIm += ftemp * phiIm;
					AmnImAcc = alphaIm + epsIm;
					epsIm += (alphaIm - AmnImAcc);

					
					/* AmnImAcc += ftemp * phiIm; */
					/* AmnRealAcc += ftemp * phiRe; */
				}
			}
			// and save the coeffs
			gsl_matrix_set(AmnReal, i, j, AmnRealAcc*dxy);
			// where is this -1 coming from? it's in the R code?
			gsl_matrix_set(AmnIm, i, j, -1.0 * AmnImAcc*dxy);
		}
	}

	gsl_matrix_free(rMat);
	gsl_matrix_free(thMat);
	gsl_vector_free(xvec);
	gsl_matrix_free(lamMat);
}

/*
 * compute the centre of mass of the grid
 */

void compute_com(gsl_matrix *array, int npts,double* cmx, double* cmy)
{
	int i, j;
	double etot = 0.0;
	double ex = 0.0, ey = 0.0;
	double xv, yv;
	/* double dr = 0.1; // size of a grid cell */
	/* double r  = 0.0; */

	// using units of -1..1 on the box
	double dx = (2*fabs(xmin))/((double)npts-1);
	

	gsl_vector *xvec = gsl_vector_alloc(npts);

	// fill in r and Theta matrices
	for(i = 0; i < npts; i++){
		gsl_vector_set(xvec ,i, xmin + dx*i);
	}
	
	for(i = 0; i < npts; i++){			
		xv = gsl_vector_get(xvec, i);
		for(j = 0; j < npts;j ++){
			yv = gsl_vector_get(xvec, j);
			ex += xv*gsl_matrix_get(array, i, j);
			ey += yv*gsl_matrix_get(array, i, j);
			etot += gsl_matrix_get(array, i, j);
		}
	}
	
	*cmx = (ex / etot);
	*cmy = (ey / etot);
}

#ifdef BUILDBIN
int main (int argc, char* argv[]){
	FILE *fptr;
	char buffer[256];
	gsl_matrix *midplane;
	gsl_matrix *AmnRe, *AmnIm;
	int i, j;
	int mmax=2, nmax=4;
	int npts = 200;
	int xtemp, ytemp;
	double etemp;

	double cmx = 0.0, cmy = 0.0;
	double dx = 2*(fabs(xmin))/((double)npts-1);

	double mod = 0.0;
	
	if(argc < 2){
		printf("# run with path to midplane file\n");
		exit(-1);
	}

	sprintf(buffer, "%s", argv[1]);
	
	printf("# reading from: %s\n", buffer);
	
	fptr = fopen(buffer, "r");
	midplane = gsl_matrix_alloc(npts, npts);
	

	while(fscanf(fptr, "%d %d %*d %lf %*f %*f %*f %*f",
							 &xtemp, &ytemp, &etemp) != EOF){
		gsl_matrix_set(midplane, xtemp - 1, ytemp -1, etemp);
	};
	
	fclose(fptr);

	AmnRe = gsl_matrix_alloc((2*mmax+1), nmax);
	AmnIm = gsl_matrix_alloc((2*mmax+1), nmax);

	compute_com(midplane, npts, &cmx, &cmy);
	
	printf("# com: %lf %lf (%d %d)\n", cmx, cmy, (int)((cmx+(fabs(xmin)))/dx), (int)((cmy+(fabs(xmin)))/dx));

	
	compute_amn(mmax, nmax, midplane, npts, AmnRe, AmnIm, 0, 0); // first compute with cm at origin of coords

	mod = 0;
	for(i = 0; i < (2*mmax+1); i++){
		for(j = 0; j < nmax; j++){
			printf("(%d %d)[%lf %lf] ", -mmax+i, j, gsl_matrix_get(AmnRe, i, j), gsl_matrix_get(AmnIm, i,j));
			mod += pow(gsl_matrix_get(AmnRe, i, j),2.0) + pow(gsl_matrix_get(AmnIm, i,j), 2.0);
		}
		printf("\n");
	}
	printf("Mod: %lf\n", mod);

	xmin = -1.0;
	
	printf("## xmin = %lf\n", xmin);
	compute_amn(mmax, nmax, midplane, npts, AmnRe, AmnIm, cmx, cmy); // first compute with cm at origin of coords
	mod = 0;
	for(i = 0; i < (2*mmax+1); i++){
		for(j = 0; j < nmax; j++){
			printf("(%d %d)[%lf %lf] ", -mmax+i, j, gsl_matrix_get(AmnRe, i, j), gsl_matrix_get(AmnIm, i,j));
			mod += pow(gsl_matrix_get(AmnRe, i, j),2.0) + pow(gsl_matrix_get(AmnIm, i,j), 2.0);
		}
		printf("\n");
	}
	printf("Mod: %lf\n", mod);
	
	gsl_matrix_free(AmnRe);
	gsl_matrix_free(AmnIm);
	gsl_matrix_free(midplane);
	return EXIT_SUCCESS;
}
#endif
