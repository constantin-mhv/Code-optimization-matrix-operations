/*
 * Tema 2 ASC
 * 2021 Spring
 */
#include "utils.h"
#include <string.h>

double* my_solver(int N, double *A, double* B) {
	printf("NEOPT SOLVER\n");
	double *temp, *left_term, temp_num;
	int i, j, k;

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
	for (i = 0; i < N; i++) {
		/** j = i because B * B' is symmetric matrix and
		 * temp[i][j] = temp [j][i], i != j */
   		for (j = i; j < N; j++){
      		for (k = 0; k < N; k++){
			/* It is no need to transpose matrix, X[i][j] = X'[j][i] */
				temp_num = B[i * N + k] * B[j * N + k];
				temp[i * N + j] += temp_num;
				/* "fill" second part of matrix */
				if (i != j)
					temp[j * N + i] += temp_num;
      		}
			
   		}
	}

	/* A * temp, result is saved in left_term */
	for (i = 0; i < N; i++){
   		for (j = 0; j < N; j++){
			/* k = i because A is upper triangular matrix */
      		for (k = i; k < N; k++)
				left_term[i * N + j] += A[i * N + k] * temp[k * N + j];
   		}
	}

	/* A' * A, result is saved in temp*/
	memset(temp, 0, N * N * sizeof(double));
	for (i = 0; i < N; i++) {
		/** j = i because A' * A is symmetric matrix and
		 * temp[i][j] = temp [j][i], i != j */
   		for (j = i; j < N; j++) {
			/* k <= i because A' is lower triangular matrix  */
			for (k = 0; k <= i; k++) {
			/* It is no need to transpose matrix, X[i][j] = X'[j][i] */
				temp_num = A[k * N + i] * A[k * N + j];
				temp[i * N + j] += temp_num;
				/* "fill" second part of matrix */
				if (i != j)
					temp[j * N + i] += temp_num;
      		}
			
   		}
	}

	/* left_term + right_term (right_term is temp), result is saved in temp */
	for (i = 0; i < N * N; i++) {
		temp[i] += left_term[i];
	}

	free(left_term);
	return temp;
}
