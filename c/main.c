#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include "cblas.h"
#include "lapacke.h"

extern int dgeqr3( int m, int n, double *A, int lda, double *T, int ldt );

int main(int argc, char ** argv) {

	int i, j, info, lda, ldq, m, n, lwork, lwork_max; 
	double *A, *Q, *As, *tau, *work=NULL;
	double normA, normR;
	double elapsed_refL, perform_refL;
	struct timeval tp;

	srand(0);

    	m = 20;
    	n = 15;
	lda = -1;
	ldq = -1;

	for(i = 1; i < argc; i++){
		if( strcmp( *(argv + i), "-ldq") == 0) {
			ldq  = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-lda") == 0) {
			lda  = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-m") == 0) {
			m  = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-n") == 0) {
			n  = atoi( *(argv + i + 1) );
			i++;
		}
	}

	if( lda < 0 ) lda = m;
	if( ldq < 0 ) ldq = m;

	printf("m = %4d, ",m);
	printf("n = %4d, ",n);
	printf("lda = %4d, ",lda);
	printf("ldq = %4d, ",ldq);
	printf("\n");

	A = (double *) malloc(lda * n * sizeof(double));
	As = (double *) malloc(lda * n * sizeof(double));
	Q = (double *) malloc(ldq * n * sizeof(double));
	tau = (double *) malloc( n * sizeof(double));

 	for(i = 0; i < lda * n; i++)
		*(A + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

 	for(i = 0; i < ldq * n; i++)
		*(Q + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, n, A, lda, As, lda );

	normA = LAPACKE_dlange_work(LAPACK_COL_MAJOR, 'F', m, n, A, lda, work );

	lwork = -1;
	work = (double *) malloc( 1 * sizeof(double));
	info = LAPACKE_dgeqrf_work( LAPACK_COL_MAJOR, m, n, A, lda, tau, work, lwork );
	printf("lwork GEQRF = %d\n",((int) work[0]));
	lwork_max = ((int) work[0]);
	lwork = -1;
	info = LAPACKE_dorgqr_work( LAPACK_COL_MAJOR, m, n, n, A, lda, tau, work, lwork );
	printf("lwork ORGQR = %d\n",((int) work[0]));
	if( ((int) work[0]) > lwork_max ) lwork_max = ((int) work[0]); 
	lwork = lwork_max;
	printf("lwork       = %d\n",lwork);

	work = (double *) malloc( lwork * sizeof(double));

	gettimeofday(&tp, NULL);
	elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

//info = LAPACKE_dgeqrf_work( LAPACK_COL_MAJOR, m, n, A, lda, tau, work, lwork );

	double *T;
	int ldt;
	ldt = n;
	T = (double *) malloc(ldt * n * sizeof(double));


	info = dgeqr3( m, n, A, lda, T, ldt );

 	for(j = 0; j < n; j++) tau[j] = T[j+j*ldt];


free(T);

	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, n, A, lda, Q, ldq );

	info = LAPACKE_dorgqr_work( LAPACK_COL_MAJOR, m, n, n, Q, ldq, tau, work, lwork );

	gettimeofday(&tp, NULL);
	elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

	free( work );
	free( tau );
	free( Q );
	free( A );

	cblas_dtrmm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, m, n, (1.0e+00), A, lda, Q, ldq );

 	for(i = 0; i < m; i++) for(j = 0; j < n; j++) As[i+j*lda] -= Q[i+j*ldq];

	normR = LAPACKE_dlange_work(LAPACK_COL_MAJOR, 'F', m, n, As, lda, work );

	perform_refL = ( 4.0e+00 * ((double) m) * ((double) n) * ((double) n) - 4.0e+00 / 3.0e+00 * ((double) n) * ((double) n) * ((double) n) )  / elapsed_refL / 1.0e+9 ;

	printf("LAPACK        :: time = %f   GFlop/sec = %f   res = %e \n", elapsed_refL, perform_refL, normR / normA );

	return 0;
}
