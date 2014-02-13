#include "fourier-decomp.h"

/**
 * setup the decomp, allocs the variable grid
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

/**
 * insert energy density val into position i,j in the matrix
 */
void fill_grid_(int *i, int *j, double *val)
{
	//printf("#(c) (%d %d) = %lf\n", *i-1, *j-1, *val);
	gsl_matrix_set(grid, *i-1, *j-1, *val);
}

/**
 * compute the fourier decomposition of a filled  grid
 * and write the results out to a file
 * 
 * this function was setup to be callable from a fortran hydro-simulator routine 
 * see also the example driver for another way to compute the coeffs
 */

void do_fdecomp_(int* mmax_in, int* nmax_in, char* outname[])
{
	int mmax;
	int nmax;
	int i,j;
	double cmx =0.0, cmy=0.0;
	FILE* fptr;
	gsl_matrix * amnRe = NULL;  
	gsl_matrix * amnIm = NULL;
  char buffer[256];

  if(*mmax_in <= 0 || *mmax_in > MAXMDECOMP){
    mmax = MMAXDEFAULT;
  } else {
    mmax = *mmax_in;
  }
  if(*nmax_in <=0 || *nmax_in > MAXNDECOMP){
    nmax = NMAXDEFAULT;
  } else {
    nmax = *nmax_in;
  }
  amnRe = gsl_matrix_alloc(2*mmax+1, nmax);
  amnIm = gsl_matrix_alloc(2*mmax+1, nmax);

	compute_com(grid, nptsGrid, &cmx, &cmy);
	
	fprintf(stderr, "# cm (%lf %lf)\n", cmx, cmy);
	fprintf(stderr, "# started computing fdecomp\n");
	compute_amn(mmax, nmax, grid, nptsGrid, amnRe, amnIm, cmx, cmy);
	
	sprintf(buffer, "%s", outname);  
	fptr = fopen(buffer, "a");
	fprintf(stderr, "# started writing fdecomp to: %s\n", buffer);
	// dump the decomp to file
	fprintf(fptr, "#\n"); // sep events with a #
	for(i = 0; i < (2*mmax+1); i++){
		for(j = 0; j < nmax; j++){
			fprintf(fptr, "%d %d %lf %lf\n", -mmax+i, j, gsl_matrix_get(amnRe, i, j), gsl_matrix_get(amnIm, i,j));
		}
	}
	fclose(fptr);
	
	gsl_matrix_free(amnRe);
	gsl_matrix_free(amnIm);
}

/**
 * Compute the fourier decomposition of array to -m:m in angular components and 1:nmax in radial components
 * 
 * Note that the gridding scheme used here is defined on [-1..1] x [-1..1], the cmx and cmy specified here should
 * be given in these units. The function "compute_com" defined below can be used to do this.
 * 
 * @arg mmax - largest angular moment computed
 * @arg nmax - largest radial moment computed
 * @arg array - 2d matrix (npts x npts) of energy density in the event
 * @arg npts - number of points in in the array
 * @arg AmnReal - 2d matrix ((2*mmax+1) x nmax), filled with Real parts of the coeffs on return
 * @arg AmnIm - 2d matrix ((2*mmax+1) x nmax), filled with Im parts of the coeffs on return
 * @arg cmx - x location of the CM of the event
 * @arg cmy - y location of the CM of the event
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

					/* kahan compensated summation (http://en.wikipedia.org/wiki/Kahan_summation_algorithm)
           * we're adding up a lot of little numbers here
           * this trick keeps accumulation errors from, well, accumulating
           */
					alphaRe = AmnRealAcc;
					epsRe += ftemp * phiRe;
					AmnRealAcc = alphaRe + epsRe;
					epsRe += (alphaRe - AmnRealAcc);

					alphaIm = AmnImAcc;
					epsIm += ftemp * phiIm;
					AmnImAcc = alphaIm + epsIm;
					epsIm += (alphaIm - AmnImAcc);
				}
			}
			// and save the coeffs
			gsl_matrix_set(AmnReal, i, j, AmnRealAcc*dxy);
			gsl_matrix_set(AmnIm, i, j, -1.0 * AmnImAcc*dxy);
		}
	}

	gsl_matrix_free(rMat);
	gsl_matrix_free(thMat);
	gsl_vector_free(xvec);
	gsl_matrix_free(lamMat);
}

/**
 * compute the centre of mass of the grid (array)
 * 
 * the grid has coordinates [-1..1] x [-1..1] and the computed CM is returned in this 2d interval
 * 
 * @arg array (npts x npts) gridded energy density of the event
 * @arg npts size of the grid
 * @arg cmx - set to the x location of the cm
 * @arg cmy - set to the y location of the cm
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



