#include "lila.h"

int main(int argc, char ** argv){

	int m, n, nb, vb, i, j, lda, ldq, ldt, info, lwork;
	double *A, *As, *T, *Q, *work=NULL;
	double normA, normR, elapsed_refL, perform_refL, norm_orth;
	struct timeval tp;

	srand(0);  // What does this do?

    	m = 30;
    	n = 20;
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
		if( strcmp( *(argv + i), "-nb") == 0) {
			nb  = atoi( *(argv + i + 1) );
			i++;
		}
	}

	if( lda < 0 ) lda = m;
	if( ldq < 0 ) ldq = m;

	printf("m = %4d; ",m);
	printf("n = %4d; ",n);
	printf("nb = %4d; ",nb);
	printf("lda = %4d; ",lda);
	printf("ldq = %4d; ",ldq);
	printf("\n");

	A = (double *) malloc(lda * n * sizeof(double));
	As = (double *) malloc(lda * n * sizeof(double));
	Q = (double *) malloc(ldq * n * sizeof(double));

 	for(i = 0; i < lda * n; i++)
		*(A + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;
 	for(i = 0; i < ldq * n; i++)
		*(Q + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, n, A, lda, As, lda );
	normA = LAPACKE_dlange_work(LAPACK_COL_MAJOR, 'F', m, n, A, lda, work );

	ldt = n;
	T = (double *) malloc(ldt * n * sizeof(double));

	lwork = 1920000;
	work = (double *) malloc( 1920000 * sizeof(double));

	gettimeofday(&tp, NULL);
	elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

	j = 0;
	if ( nb > n ) vb = n; else vb = nb;
	info = lila_dgeqr2_w02a( m, vb, 0, ldt, A, lda, T, ldt, Q, ldq, work, lwork );
	j += vb;
	if ( j+nb > n ) vb = n-j; else vb = nb;
	while( vb!=0 ){
	info = lila_dormqrf_w02( m, vb, j, 0, j, ldt, A, lda, T, ldt, work, lwork );
	info = lila_dgeqr2_w02a( m, vb, j, ldt, A, lda, T, ldt, Q, ldq, work, lwork );
	info = lila_dlarft_connect_w02(m, vb, j, 0, -1, A, lda, T, ldt );
	info = lila_dormqrbz_w02( m, vb, j, 0, j, ldt, A, lda, Q, ldq, T, ldt, work, lwork );	
	j += vb;
	if ( j+nb > n ) vb = n-j; else vb = nb;
	}

	gettimeofday(&tp, NULL);
	elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

	free( work );
	free( T );

	// This check is for the orthogonality of Q
	work = (double *) malloc(n * n * sizeof(double));
	info = LAPACKE_dlaset( LAPACK_COL_MAJOR, 'A', n, n, (0e+00), (1e+00), work, n );
	cblas_dsyrk( CblasColMajor, CblasUpper, CblasTrans, n, m, 1.0e+00, Q, ldq, -1.0e+00, work, n );
	norm_orth = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', n, n, work, n, NULL );
	free( work );

	// This check to verify A=QR
	cblas_dtrmm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, m, n, (1.0e+00), A, lda, Q, ldq );
 	for(i = 0; i < m; i++) for(j = 0; j < n; j++) As[i+j*lda] -= Q[i+j*ldq];
	normR = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', m, n, As, lda, work );
	perform_refL = ( 4.0e+00 * ((double) m) * ((double) n) * ((double) n) - 4.0e+00 / 3.0e+00 * ((double) n) * ((double) n) * ((double) n) )  / elapsed_refL / 1.0e+9 ;
	printf("| time = %f   GFlop/sec = %f   res = %e    orth = %e \n", elapsed_refL, perform_refL, normR / normA, norm_orth );

	free( Q );
	free( A );
	free( As );

	return 0;
}
