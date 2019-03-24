#include "lila.h"

int main(int argc, char ** argv) {

	double *tmp;
	int lda, ldq, m, n, i, j, info, lwork;
	double normA, normR;
	double *A, *As, *tau, *Q, *work;
	double elapsed_refL;
	struct timeval tp;
	double norm_orth;

	srand(0);

	m = 50;
	n = 10;
	lda = 51;
	ldq = 52;

	for(i = 1; i < argc; i++){
		if( strcmp( *(argv + i), "-m") == 0) {
			m  = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-n") == 0) {
			n  = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-lda") == 0) {
			lda  = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-ldq") == 0) {
			ldq  = atoi( *(argv + i + 1) );
			i++;
		}
	}

	if (lda < m) lda=m;
	if (ldq < m) ldq=m;
	
	//printf("m = %4d; ",m);
	//printf("n = %4d; ",n);
	//printf("lda = %4d; ",lda);
	
	printf(" %5d %5d ",m,n);

	A  = (double *) malloc(lda * n * sizeof(double));
	As  = (double *) malloc(lda * n * sizeof(double));

	tau  = (double *) malloc(n * sizeof(double));
	Q  = (double *) malloc(ldq * n * sizeof(double));

 	for(i = 0, tmp=A; i < lda * n; i++, tmp++){
		*tmp = (double)rand() / (double)(RAND_MAX) - 0.5e+00;
	}

	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, n, A, lda, As, lda );
	normA = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', m, n, A, lda, NULL );	


	work = (double *) malloc(1 * sizeof(double));
	lwork = -1;
	info = LAPACKE_dgeqrf_work( LAPACK_COL_MAJOR, m, n, A, lda, tau, work, -1);
//	printf("work[0] = %f;", work[0] );
	lwork = ((int) work[0]);
	info = LAPACKE_dorgqr_work( LAPACK_COL_MAJOR, m, n, n, Q, ldq, tau, work, -1);
//	printf(" work[0] = %f;", work[0] );
	if( lwork < ((int) work[0]) ) lwork = ((int) work[0]);
	free( work );
	work = (double *) malloc(lwork * sizeof(double));


	gettimeofday(&tp, NULL);
	elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
	info = LAPACKE_dgeqrf_work( LAPACK_COL_MAJOR, m, n, A, lda, tau, work, lwork);
	LAPACKE_dlacpy( LAPACK_COL_MAJOR, 'A', m, n, A, lda, Q, ldq );
	info = LAPACKE_dorgqr_work( LAPACK_COL_MAJOR, m, n, n, Q, ldq, tau, work, lwork);
	gettimeofday(&tp, NULL);
	elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
//	printf(" time = %f;", elapsed_refL );
//	printf(" perf = %f;", ( 4.0e+00 * ((double) m) * ((double) n) * ((double) n) - 4.0e+00 / 3.0e+00 * ((double) n) * ((double) n) * ((double) n) )  / elapsed_refL / 1.0e+9 );
	printf(" %f", elapsed_refL );
	printf(" %f\n", ( 4.0e+00 * ((double) m) * ((double) n) * ((double) n) - 4.0e+00 / 3.0e+00 * ((double) n) * ((double) n) * ((double) n) )  / elapsed_refL / 1.0e+9 );

	free(work);
	free(tau);

	return 0;



	work = (double *) malloc(n * n * sizeof(double));

	info = LAPACKE_dlaset( LAPACK_COL_MAJOR, 'A', n, n, (0e+00), (1e+00), work, n );

	cblas_dsyrk( CblasColMajor, CblasUpper, CblasTrans, n, m, 1.0e+00, Q, ldq, -1.0e+00, work, n );

	norm_orth = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', n, n, work, n, NULL );

	cblas_dtrmm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, m, n, (1.0e+00), A, lda, Q, ldq );

 	for(i = 0; i < m; i++) for(j = 0; j < n; j++) As[i+j*lda] -= Q[i+j*ldq];

	normR = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', m, n, As, lda, work );

	printf(" res = %e    orth = %e \n", (normR / normA), norm_orth );

	free(Q);
	free(A);
	free(As);
	free( work );

	return 0;
}
