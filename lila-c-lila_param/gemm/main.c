#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>

#include "cblas.h"

int main(int argc, char ** argv) {

	int i, info, m, n, k;
	double *A, *B, *C;
	double elapsed_ref, perform_ref;
	struct timeval tp;

	srand(0);

    	m = 1000;
    	n = m;
	k = m;

	for(i = 1; i < argc; i++){
		if( strcmp( *(argv + i), "-m") == 0) {
			m  = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-n") == 0) {
			n  = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-k") == 0) {
			k  = atoi( *(argv + i + 1) );
			i++;
		}
	}

	A = (double *) malloc( m * k * sizeof(double));
	B = (double *) malloc( k * n * sizeof(double));
	C = (double *) malloc( m * n * sizeof(double));

 	for(i = 0; i < m * k; i++)
		*(A + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

 	for(i = 0; i < k * n; i++)
		*(B + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

 	for(i = 0; i < m * n; i++)
		*(C + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

	// Begin timing
	gettimeofday(&tp, NULL);
	elapsed_ref=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

	cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, k, (+1.0e+00), A, m, B, n, (+1.0e+00), C, n);	

	gettimeofday(&tp, NULL);
	elapsed_ref+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

	perform_ref = ( 2.0e+00 * ((double) m) * ((double) k) * ((double) n) ) / elapsed_ref / 1.0e+9 ;

	printf("dgemm - ");
	printf("m = %4d, ",         m);
	printf("n = %4d, ",         n);
	printf(" time = %f    GFlop/sec = %f ", elapsed_ref, perform_ref);	
	printf(" \n");

	free( A );
	free( B );
	free( C );

	return 0;
}
