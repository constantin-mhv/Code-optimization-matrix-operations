/*
 * Tema 2 ASC
 * 2021 Spring
 */
#include "utils.h"
#include <cblas.h>
#include <string.h>

double* my_solver(int N, double *A, double *B) {
	printf("BLAS SOLVER\n");
	double *temp, *left_term;

	temp = malloc(N * N * sizeof(double));
	if (!temp) {
		fprintf(stderr, "Memory allocation ERROR\n");
		return NULL;
	}

	left_term = malloc(N * N * sizeof(double));
	if (!left_term) {
		fprintf(stderr, "Memory allocation ERROR\n");
		free(temp);
		return NULL;
	}

	/* C = A * B * B' + A' * A */
	
	/* A * B, result is saved in temp */
	memcpy(temp, B, N * N * sizeof(double));
	cblas_dtrmm(
		CblasRowMajor,
		CblasLeft,
		CblasUpper,
		CblasNoTrans,
		CblasNonUnit,
		N, N,
		1.0, A, N,
		temp, N
	);

	/* temp * B', result is saved in left_term*/
	cblas_dgemm(
		CblasRowMajor,
		CblasNoTrans,
		CblasTrans,
		N, N, N,
		1.0, temp, N,
		B, N, 1.0,
		left_term, N
	);

	/* A' * A, result is saved in temp. Now temp is right term */
	memcpy(temp, A, N * N * sizeof(double));
	cblas_dtrmm(
		CblasRowMajor,
		CblasLeft,
		CblasUpper,
		CblasTrans,
		CblasNonUnit,
		N, N,
		1.0, A, N,
		temp, N
	);

	/* left_term + right_term (right_term is temp), result is saved in temp */
	cblas_daxpy(
		N * N,
		1.0, left_term, 1,
		temp, 1
	);

	free(left_term);
	return temp;
}
