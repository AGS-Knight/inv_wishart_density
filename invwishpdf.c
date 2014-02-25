#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_permutation.h>
#include "invwishpdf.h"

int main(int argc, char *argv[]){

	int args;
	int m, n;
	double pdf, dof;

	args = 1;

	m = atoi(argv[args++]);
	n = atoi(argv[args++]);
	
	// allocate space for matrix X, test_copy, and Scale.
	gsl_matrix *X = gsl_matrix_alloc(m, n);
	gsl_matrix *inv = gsl_matrix_alloc(m,n);
	gsl_matrix *Scale = gsl_matrix_alloc(m, n);

	gsl_matrix_set_identity(Scale);
	fill_matrix(X);

	dof = m + 1;
	pdf = iwishpdf(X, Scale, inv, dof);

	// Print the data
	printf("\n%s", "Inverse Wishart matrix");
	print_matrix(X);
	printf("\n%s", "Scale matrix");
	print_matrix(Scale);
	printf("\n%s", "Wishart matrix");
	print_matrix(inv);
	printf("\n%s", "Density");
	printf("\n%f\n", pdf);

	return(0);
}

gsl_matrix fill_matrix(gsl_matrix *X)
{
	int i,j;
	int step = 0;
	int n = X->size1;
	int m = X->size2;

	// Create a test matrix
	double a_vect[] = {0.17955079, -0.09332222, -0.01273225, -0.09332222,  0.28282278,  0.01724885, -0.01273225,  0.01724885,  0.14424177};
	for (i = 0; i < m; ++i){
		for (j = 0; j < n; ++j){
			gsl_matrix_set(X, i, j, a_vect[step++]);
		}
	}
	return(*X);	
}

double iwishpdf(gsl_matrix *X, gsl_matrix *Scale, gsl_matrix *inv, double dof)
{
	double X_det, scale_det, denom, pdf, trace, numerator;	
	int m = X->size1;
	int n = X->size2;
	
	// Allocate matrix for inv(X)
	gsl_matrix *for_mult = gsl_matrix_alloc(m, n);
	
	// Get determinant of X and Scale matrix
	X_det = matrix_determ(X);
	scale_det = matrix_determ(Scale);
	
	// Invert X
	inv_matrix(X,inv);
	
	// Multiple Scale * inv(X)
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Scale, inv, 1.0, for_mult);
	
	// Get trace of above.
	trace = matrix_trace(for_mult);

	numerator =  pow(scale_det, dof / 2.0) * pow(X_det, (-dof-m-1)/ 2.0) * exp(-0.5 * trace);
	denom = pow(2,dof * m / 2) * mv_gamma(dof/2, m);

	pdf = (numerator/denom);

	return(pdf);
}

gsl_matrix inv_matrix(gsl_matrix *X, gsl_matrix *inv)
{
	int s;
	int n = X->size1;
	int m = X->size2;
	gsl_matrix *a_copy = gsl_matrix_alloc(m, n);
	gsl_matrix_memcpy( a_copy, X );
	gsl_permutation *P = gsl_permutation_alloc(n);
	gsl_linalg_LU_decomp(a_copy, P, &s);
	gsl_linalg_LU_invert(a_copy, P, inv);
	return(*inv);
}


double matrix_determ(gsl_matrix *X)
{
	int s;
	int n = X->size1;
	int m = X->size2;
	gsl_matrix *a_copy = gsl_matrix_alloc(m, n);
	gsl_matrix_memcpy(a_copy, X );
	gsl_permutation *P = gsl_permutation_alloc(n);
	gsl_linalg_LU_decomp(a_copy, P, &s);
	double my_det = gsl_linalg_LU_det (a_copy, s);
	return(my_det);
}

double matrix_trace(gsl_matrix *X){
	int i, m;
	m = X->size1;
	double trace = 0.0;
	for(i=0;i<m;i++){
		trace += gsl_matrix_get(X,i,i);
	}
	return(trace);
}

double mv_gamma(double a, double d){
	double val = 1.0;
	int i;
	for(i = 1; i <= d; i++){
		val *= gsl_sf_gamma(a - (0.5 * (i - 1)));
	}
	val *=  pow(M_PI, (d * (d - 1) / 4.0));
	return(val);
}

void print_matrix(gsl_matrix *X)
{
	int i, j;
	int n = X->size1;
	int m = X->size2;
	printf("\n");
	for(i=0; i<m; i++){
		for(j=0; j<n; j++){
			printf("%05.2f ", gsl_matrix_get(X, i, j));
		}
		printf("\n");
	}
}
