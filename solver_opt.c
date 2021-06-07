/*
 * Tema 2 ASC
 * 2021 Spring
 */
#include "utils.h"
#include <string.h>

double* my_solver(int N, double *A, double* B) {
	printf("OPT SOLVER\n");
	double *temp, *left_term;
	/** temp_i_j and left_term_i_j will be used to write calculated result
	 * more efficient
	 * orig_pa, pa, pb will be used to access more efficient data from matrices
	*/
	double *orig_pa, *pa, *pb, *temp_i_j, *left_term_i_j;
	int i, j, k;
	register double sum = 0;

	temp = calloc(N * N, sizeof(double));
	if (! temp) {
		fprintf(stderr, "Memory allocation ERROR\n");
		return NULL;
	}

	left_term = calloc(N * N, sizeof(double));
	if (!left_term) {
		fprintf(stderr, "Memory allocation ERROR\n");
		free(temp);
		return NULL;
	}

	/* C = A * B * B' + A' * A */
	
	/* B * B', result is saved in temp */
	temp_i_j = &temp[0];
	for (i = 0; i < N; ++i, temp_i_j += N){
		orig_pa = &B[i * N];
		/** j = i because B * B' is symmetric matrix and
		 * temp[i][j] = temp [j][i], i != j */
   		for (j = i; j < N; ++j) {
			pa = orig_pa;
			pb = &B[j * N];
			sum = 0;
      		for (k = 0; k < N; ++k){
			/* It is no need to transpose matrix, X[i][j] = X'[j][i] */
				sum += *pa * *pb;
				pa++;
				pb ++;
      		}
			*(temp_i_j + j) = sum;
			/* "fill" second part of matrix */
			if (i != j)
				temp[j * N + i] = sum;
   		}
	}

	/* A * temp, result is saved in left_term */
	left_term_i_j = &left_term[0];
	for (i = 0; i < N; ++i, left_term_i_j += N){
		orig_pa = &A[i * N];
   		for (j = 0; j < N; ++j) {
			pa = orig_pa + i;
			pb = &temp[j] + N * i;
			sum = 0;
			/* k = i because A is upper triangular matrix */
      		for (k = i; k < N; ++k){
				sum += *pa * *pb;
				pa++;
				pb+=  N;
      		}
			*(left_term_i_j + j) = sum;
			// left_term[i * N + j] = sum;
   		}
	}

	/* A' * A, result is saved in temp*/
	memset(temp, 0, N * N * sizeof(double));
	temp_i_j = &temp[0];
	for (i = 0; i < N; ++i, temp_i_j += N) {
		/** j = i because A' * A is symmetric matrix and
		 * temp[i][j] = temp [j][i], i != j */
   		for (j = i; j < N; ++j) {
			pa = &A[i];
			pb = &A[j];
			sum = 0;
			/* k <= i because A' is lower triangular matrix  */
			for (k = 0; k <= i; ++k) {
			/* It is no need to transpose matrix, X[i][j] = X'[j][i] */
				sum += *pa * *pb;
				pa += N;
				pb += N;
      		}
			*(temp_i_j + j) = sum;
			/* "fill" second part of matrix */
			if (i != j)
			 	temp[j * N + i] = sum;
			
   		}
	}

	/* left_term + right_term (right_term is temp), result is saved in temp */
	pa = &temp[0];
	pb = &left_term[0];
	for (i = 0; i < N * N; ++i) {
		*pa += *pb;
		pa++;
		pb++;
	}
	
	free(left_term);
	return temp;
}
