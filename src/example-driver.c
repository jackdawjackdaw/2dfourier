#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "fourier-decomp.h"

/**
 * ccs, 13.02.2014
 * example driver routine, reads an input file and computes the decomposition moments
 * Amn
 */ 

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
