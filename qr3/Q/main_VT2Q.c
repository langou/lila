#include "qr3.h"

int main(int argc, char ** argv) {

	int i, info, lda, ldt, ldq, m, n, verbose, testing;
	int lwork;
	double *A, *As, *T, *Q, *tau, *work=NULL;
	double orth, repres_1;
	double elapsed_ref1, perform_ref1, elapsed_ref2, perform_ref2;
	struct timeval tp;
	
	srand(0);

    	m         = 27;
    	n         = 17;
	lda       = -1;
	ldq       = -1;
	verbose   = 0;
	testing   = 1;

	for(i = 1; i < argc; i++){
		if( strcmp( *(argv + i), "-lda") == 0) {
			lda  = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-ldq") == 0) {
			ldq  = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-verbose") == 0) {
			verbose  = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-testing") == 0) {
			testing  = atoi( *(argv + i + 1) );
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

	if( m < n ){ printf("\n\n We only work with m >=n, please reconsider your m and n\n\n"); return 0; }

	if( lda < 0 ) lda = m;
	if( ldq < 0 ) ldq = m;

	ldt = n;

	A  = (double *) malloc( lda * n * sizeof(double));
 	As = (double *) malloc( lda * n * sizeof(double));
	T  = (double *) malloc( ldt * n * sizeof(double));
	Q  = (double *) malloc( ldq * n * sizeof(double));

 	for(i = 0; i < lda * n; i++)
		*(A + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

	for(i = 0; i < ldq * n; i++)
		*(Q + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, n, A, lda, As, lda );

	tau   = (double *) malloc( n * sizeof(double));
	lwork = -1;
	work  = (double *) malloc( 1 * sizeof(double));
	info = LAPACKE_dgeqrf_work( LAPACK_COL_MAJOR, m, n, A, lda, tau, work, lwork ); 
	lwork = (int) work[0];
	free(work);
	work  = (double *) malloc( lwork * sizeof(double));

	gettimeofday(&tp, NULL);
	elapsed_ref1=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
	info = LAPACKE_dgeqrf_work( LAPACK_COL_MAJOR, m, n, A, lda, tau, work, lwork ); 
	gettimeofday(&tp, NULL);
	elapsed_ref1+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
	perform_ref1 = ( 2.0e+00 * ((double) m) * ((double) n) * ((double) n) - 2.0e+00 / 3.0e+00 * ((double) n) * ((double) n) * ((double) n) )  / elapsed_ref1 / 1.0e+9 ;

	free( work );

	info = LAPACKE_dlarft_work( LAPACK_COL_MAJOR, 'F', 'C', m, n, A, lda, tau, T, ldt);

	gettimeofday(&tp, NULL);
	elapsed_ref2=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);


	// LAPACK version for constructing Q
	//lwork = -1;
	//work  = (double *) malloc( 1 * sizeof(double));
	//info = LAPACKE_dorgqr_work( LAPACK_COL_MAJOR, m, n, n, Q, ldq, tau, work, lwork );
	//lwork = (int) work[0];
	//free(work);
	//work  = (double *) malloc( lwork * sizeof(double));

	//info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, n, A, lda, Q, ldq );
	//info = LAPACKE_dlaset( LAPACK_COL_MAJOR, 'U', n, n, (0.0e+00), (1.0e+00), Q, ldq);
	//info = LAPACKE_dorgqr_work( LAPACK_COL_MAJOR, m, n, n, Q, ldq, tau, work, lwork );
	//free(work);

	// Setting ourselves up.
	// the lower part of Q is V, and the upper part of Q is T
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'L', m-1, n, A+1, lda, Q+1, ldq );
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'U', n, n, T, ldt, Q, ldq );

	info = ULTinU( n, Q, ldq, Q, ldq );

	cblas_dtrmm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, m-n, n, (-1.0e+00), Q, ldq, Q+n, ldq );

	info = mLUinA( n, Q, ldq );

	for(i = 0; i < n; i++) Q[ i + ldq * i ] = 1.00e+00 + Q[ i + ldq * i ];

	gettimeofday(&tp, NULL);
	elapsed_ref2+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
	perform_ref2 = ( 2.0e+00 * ((double) m) * ((double) n) * ((double) n) - 2.0e+00 / 3.0e+00 * ((double) n) * ((double) n) * ((double) n) )  / elapsed_ref2 / 1.0e+9 ;

	free( tau );

	if ( verbose ){ 

		printf("ORGQR_""recursive"" - ");
		printf("m = %4d, ",         m);
		printf("n = %4d, ",         n);
		printf("lda = %4d, ",     lda);
		printf(" \n");
		printf(" time = %f    GFlop/sec = %f ", elapsed_ref1, perform_ref1);	
		printf(" timeT = %f    GFlop/secT = %f ", elapsed_ref2, perform_ref2);	
		printf(" \n ");

	} else {

		printf("%6d %6d %6d %16.8f %10.3f %16.8f %10.3f ", m, n, lda, elapsed_ref1, perform_ref1, elapsed_ref2, perform_ref2);

	} 

	if ( testing ){

		info = qr3_test_qq_orth_1( &orth, m, n, Q, ldq );		
		if ( verbose ) printf("qq_orth  = %5.1e  \n ",orth); else printf(" %5.1e  ",orth); 

		info = qr3_test_qr_repres_1( &repres_1, m, n, As, lda, Q, ldq, A, lda );
		if ( verbose ) printf("qr_repres = %5.1e  \n ",repres_1); else printf(" %5.1e  ",repres_1); 

	}

	if ( !verbose ) printf("\n");		

	free( A );
	free( As );
	free( Q );
	free( T );

	return 0;

}