#include "qr3.h"

extern void our_dorgqr_fortran_( int *m, int *n, int *k, double *a, int *lda, double *tau, double *work, int *lwork, int *info );

int main(int argc, char ** argv) {

	int i, lda, ldq, ldr, ldt, m, n, k, nb, verbose, testing;
	int lwork;
	double *A, *Q, *R, *T, *tau, *work;
	double orth, repres;
	double elapsed_ref, perform_ref;
	struct timeval tp;
	
	srand(0);

    	m         = 27;
    	k         = 17;
    	n         = 20;
	lda       = -1;
	ldq       = -1;
	ldr       = -1;
	ldt       = -1;
	nb        = 10;
	verbose   = 0;
	testing   = 1;

	for(i = 1; i < argc; i++){
		if( strcmp( *(argv + i), "-lda") == 0) {
			lda = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-ldq") == 0) {
			ldq = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-ldr") == 0) {
			ldr = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-nb") == 0) {
			nb = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-verbose") == 0) {
			verbose = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-testing") == 0) {
			testing = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-m") == 0) {
			m = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-n") == 0) {
			n = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-k") == 0) {
			k = atoi( *(argv + i + 1) );
			i++;
		}
	}

	if (( m < n )||( n < k )) { printf("\n\n We need k <= n <= m\n\n"); return 0; }

	if( lda < 0 ) lda = m;
	if( ldq < 0 ) ldq = m;
	if( ldr < 0 ) ldr = k;
	if( ldt < 0 ) ldt = k;

	A = (double *) malloc( lda * k * sizeof(double));
	Q = (double *) malloc( ldq * n * sizeof(double));
 	R = (double *) malloc( ldr * k * sizeof(double));
 	T = (double *) malloc( ldt * k * sizeof(double));

 	for(i = 0; i < lda * k; i++)
		*(A + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

	for(i = 0; i < ldq * n; i++)
		*(Q + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

	for(i = 0; i < ldr * k; i++)
		*(R + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

	tau = (double *) malloc( k * sizeof(double));

	work = (double *) malloc( 1 * sizeof(double));
	LAPACKE_dgeqrf_work( LAPACK_COL_MAJOR, m, k, Q, ldq, tau, work, -1 ); 
	lwork = ((int) work[0]);
	LAPACKE_dorgqr_work( LAPACK_COL_MAJOR, m, n, k, Q, ldq, tau, work, -1 );
	if (lwork < ((int) work[0])) lwork = ((int) work[0]); 
	free( work );
	work = (double *) malloc( lwork * sizeof(double));

	LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, k, A, lda, Q, ldq );

	gettimeofday(&tp, NULL);
	elapsed_ref=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

//	LAPACKE_dgeqrf_work( LAPACK_COL_MAJOR, m, k, Q, ldq, tau, work, lwork ); 
//	our_dgeqrf( m, k, nb, Q, ldq, tau, work, lwork );
//	dgeqr3R( m, k, Q, ldq, T, ldt, R, ldr );
//	for(i=0;i<k;i++){tau[i] = T[i+i*ldt];}
	dgeqr3R( m, k, Q, ldq, Q, ldq, R, ldr );
	for(i=0;i<k;i++){tau[i] = Q[i+i*ldq];}
//	dgeqr3R( m, k, Q, ldq, T, ldt, Q, ldq );
//	for(i=0;i<k;i++){tau[i] = T[i+i*ldt];}
//	LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'U', k, k, Q, ldq, R, ldr );

	//qr3_larft( m, k, Q, ldq, T, ldt, tau );
//	LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'U', k, k, Q, ldq, R, ldr );
//	LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'L', m-1, k, A+1, lda, Q+1, ldq );
//	LAPACKE_dlaset( LAPACK_COL_MAJOR, 'U', k, k, (3.0e+00), (2.0e+00), Q, ldq);

//	LAPACKE_dorgqr_work( LAPACK_COL_MAJOR, m, n, k, Q, ldq, tau, work, lwork );
//	int info; our_dorgqr_fortran_( &m, &n, &k, Q, &ldq, tau, work, &lwork, &info );
	//our_dorgqr( m, n, k, nb, Q, ldq, tau, work, lwork );
//	qr3_dorgqr( m, k, Q, ldq, T, ldt, tau );
//	qr3_dorgqr( m, k, Q, ldq, Q, ldq, tau );
	qr3_dorgqr_level1( m, n, k, nb, Q, ldq, tau, work, lwork );
//	qr3_dorgqr_level1_UT( m, n, k, nb, Q, ldq, tau, work, lwork );

	gettimeofday(&tp, NULL);
	elapsed_ref+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

//	dV2Q( m, n, n, nb, nb, 0, Q, ldq, tau, work, lwork, info );

	free( tau );
	free( work );

	perform_ref = ((double) flops_lapack_org2r( m, n, k ) + (double) flops_lapack_geqr2( m, k ) ) / elapsed_ref / 1.0e+9 ;


	if ( verbose ){ 

		printf("ORGQR");
		printf("m = %4d, ",         m);
		printf("n = %4d, ",         n);
		printf("k = %4d, ",         k);
		printf("nb = %4d, ",       nb);
		printf(" \n");
		printf(" time = %f    GFlop/sec = %f ", elapsed_ref, perform_ref);	
		printf(" \n ");

	} else {

		printf("%6d %6d %6d %6d %16.8f %10.3f ", m, n, k, nb, elapsed_ref, perform_ref);

	} 

	if ( testing ){

		qr3_test_qq_orth_1( &orth, m, n, Q, ldq );		
		if ( verbose ) printf("qq_orth  = %5.1e  \n ",orth); else printf(" %5.1e  ",orth); 

		qr3_test_qr_repres_1( &repres, m, k, A, lda, Q, ldq, R, ldr );
		if ( verbose ) printf("qr_repres = %5.1e  \n ",repres); else printf(" %5.1e  ",repres); 

	}

	if ( !verbose ) printf("\n");		

	free( A );
	free( Q );
	free( R );
	free( T );

	return 0;

}
