#include "qr3.h"

// CHECK PERFORM_REF --- USING LAPACK AND HAVE NOT LOOKED INTO IT YET 

int main(int argc, char ** argv) {

	int i, lda, ldq, ldr, ldt, lwork, m, n, k, nb, verbose, testing;
	double *A, *Q, *R, *T, *tau, *work;
	double orth, repres;
	double elapsed_ref, perform_ref;
	struct timeval tp;
	
	srand(0);

    	m         = 27;
    	n         = 20;
    	k         = 10;
    	nb        = 20;
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
		if( strcmp( *(argv + i), "-nb") == 0) {
			nb = atoi( *(argv + i + 1) );
			i++;
		}
	}

	if ( m < n ) { printf("\n\n We need n <= m\n\n"); return 0; }

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

	lwork = n;
	work = (double *) malloc( lwork * n * sizeof(double));
	tau  = (double *) malloc( k * sizeof(double));
	LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, k, A, lda, Q, ldq );

	dgeqr3    ( m, k, Q, ldq, T, ldt );
//	dgeqr3_ISW( m, k, Q, ldq, T, ldt ); 
//	dgeqr3R   ( m, k, Q, ldq, T, ldt, R, ldr ); // Needs nothing done ( but the allocation for space )
//	dgeqr3R   ( m, k, Q, ldq, T, ldt, Q, ldq ); // Needs copy 'U' of Q in R && diag( T ) into tau  ----  same as qr3, qr3_ISW
//	dgeqr3R   ( m, k, Q, ldq, Q, ldq, R, ldr ); // Needs to copy diag( Q ) into tau

	LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'U', k, k, Q, ldq, R, ldr );
	LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'U', k, k, T, ldt, Q, ldq );
	for(i=0;i<k;i++) tau[i] = T[i+i*ldt];
//	for(i=0;i<k;i++) tau[i] = Q[i+i*ldq];

	gettimeofday(&tp, NULL);
	elapsed_ref=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);


	qr3_null_dorgqr( m, n, k, nb, Q, ldq, Q+k*ldq, ldq, tau, work, lwork );
//	dorgqr_after( m, n, k, Q, ldq, Q, ldq, Q+k*ldq, ldq );


	gettimeofday(&tp, NULL);
	elapsed_ref+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

	qr3_dorgqr( m, k, Q, ldq, T, ldt, tau );

	free( tau );
	free( work );

//	THIS IS NOT CORRECT AT ALL >> perform_ref SPIT OUT A NEGATIVE VALUE FOR  -m 1098 -n 1098 -nb 147 -k 698   ----    nb is not included
//	perform_ref = ((double) flops_lapack_org2r( m, n, k ) + (double) flops_lapack_geqr2( m, k ) ) / elapsed_ref / 1.0e+9 ;
	long int flops_lapack_geqr2;
	long int flops_lapack_org2r;
	flops_lapack_org2r = (( 6*(m-k)*k + 4*k*k - 3*(m-k) - 1 )*k / 3 ) + ( 4*m - 2*k + 1 )*k*( n - k );
	flops_lapack_geqr2 = (( 6*m*k*k - 2*k*k*k + 3*m*k + 17*k ) / 3 );
	perform_ref = ( ((double) flops_lapack_org2r ) + ((double) flops_lapack_geqr2 ) ) / elapsed_ref / 1.0e+9 ;

	if ( verbose ){ 

		printf("ORGQR");
		printf("m = %4d, ",           m);
		printf("n = %4d, ",           n);
		printf("k = %4d, ",           k);
		printf("nb = %4d, ",         nb);
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
