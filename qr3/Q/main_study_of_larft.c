#include "qr3.h"

int main(int argc, char ** argv) {

	int i, lda, ldq, ldr, ldt, m, n, k, verbose, testing;
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
		if( strcmp( *(argv + i), "-ldt") == 0) {
			ldt = atoi( *(argv + i + 1) );
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

	k = n;

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

	for(i = 0; i < ldt * k; i++)
		*(R + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

	tau = (double *) malloc( k * sizeof(double));

	work = (double *) malloc( 1 * sizeof(double));
	LAPACKE_dgeqrf_work( LAPACK_COL_MAJOR, m, k, Q, ldq, tau, work, -1 ); 
	lwork = ((int) work[0]);
	LAPACKE_dorgqr_work( LAPACK_COL_MAJOR, m, n, k, Q, ldq, tau, work, -1 );
	if (lwork < ((int) work[0])) lwork = ((int) work[0]); 
	free( work );
	work = (double *) malloc( lwork * sizeof(double));

	gettimeofday(&tp, NULL);
	elapsed_ref=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

	LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, k, A, lda, Q, ldq ); //do we want to time this copy?

	LAPACKE_dgeqrf_work( LAPACK_COL_MAJOR, m, k, Q, ldq, tau, work, lwork ); 
//	dgeqr3_right( m, k, Q, ldq, T, ldt );
//	dgeqr3_left( m, k, Q, ldq, T, ldt );

	LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'U', k, k, Q, ldq, R, ldr );

	gettimeofday(&tp, NULL);
	elapsed_ref=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

//	LAPACKE_dlarft_work( LAPACK_COL_MAJOR, 'F', 'C', m, k, Q, ldq, tau, T, ldt );

//	qr3_larft( m, k, Q, ldq, T, ldt, tau );

//	{ int i,j; for(i=0;i<k;i++){ for(j=0;j<i;j++){ T[j+i*ldt] = Q[i+j*ldq];}} }
//	dV2N( k, T, ldt );
//	cblas_dsyrk( CblasColMajor, CblasUpper, CblasTrans, k, m-k, (+1.0e+00), Q+k, ldq, (+1.0e+00), T, ldt );
//	dN2T( k, tau, T, ldt );

	{ int i,j; for(i=0;i<k;i++){ for(j=0;j<i;j++){ T[j+i*ldt] = Q[i+j*ldq];}} }
	dV2N( k, T, ldt );
	cblas_dsyrk( CblasColMajor, CblasUpper, CblasTrans, k, m-k, (+1.0e+00), Q+k, ldq, (+1.0e+00), T, ldt );
	{ int i; for(i=0;i<k;i++){ T[i+i*ldt] = (+1.0e+00) / tau[i];}} 
	LAPACKE_dtrtri( LAPACK_COL_MAJOR, 'U', 'N', k, T, ldt );

	gettimeofday(&tp, NULL);
	elapsed_ref+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

	for(i=0;i<k;i++){tau[i] = T[i+i*ldt];}
	qr3_dorgqr( m, k, Q, ldq, T, ldt, tau );

//	LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'U', k, k, T, ldt, Q, ldq );
//	dVT2Q( m, k, Q, ldq );


	free( tau );
	free( work );

//	perform_ref = ((double) flops_larft( m, n )) / elapsed_ref / 1.0e+9 ;
	long int flops_larft;
	flops_larft = (( n*( n - 1 )*( 6*m - ( 2*n - 1 ) )  ) /  6);
	perform_ref = ((double) flops_larft) / elapsed_ref / 1.0e+9 ;
	

	if ( verbose ){ 

		printf("ORGQR");
		printf("m = %4d, ",         m);
		printf("n = %4d, ",         n);
		printf("k = %4d, ",         k);
		printf(" \n");
		printf(" time = %f    GFlop/sec = %f ", elapsed_ref, perform_ref);	
		printf(" \n ");

	} else {

		printf("%6d %6d %6d %16.8f %10.3f ", m, n, k, elapsed_ref, perform_ref);

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
	free( T );
	free( R );

	return 0;

}
