#include "../src/qr2.h"
#include "../check/check.h"
#include "../flops/flops.h"

int main(int argc, char ** argv) {

	int i, lda, ldq, ldr, ldt, m, n, k, nb, verbose, testing;
	int lwork;
	double *A, *Q, *R, *T, *tau, *work;
	double orth, repres;
	double elapsed, perform;
	struct timeval tp;
	
	srand(0);

    	m         = 27;
    	k         = 17;
    	n         = 20;
    	nb        =  3;
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
	work = (double *) malloc( lwork * sizeof(double));

	LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, k, A, lda, Q, ldq );

//	LAPACKE_dgeqrf_work( LAPACK_COL_MAJOR, m, k, Q, ldq, tau, work, lwork ); 

	lapack_ref_dgeqrf( m, k, nb, Q, ldq, tau, work, lwork );

	LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'U', k, k, Q, ldq, R, ldr );

	free( work );
	work = (double *) malloc( n * nb * sizeof(double));

	gettimeofday(&tp, NULL);

	elapsed=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

	lapack_ref_dorgqr( m, n, k, nb, Q, ldq, tau, work, lwork );

	gettimeofday(&tp, NULL);
	elapsed+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

	free( tau );
	free( work );

	perform = ((double) flops_lapack_org2r( m, n, k ) + (double) flops_lapack_geqr2( m, k ) ) / elapsed / 1.0e+9 ;

	if ( verbose ){ 

		printf("ORGQR  ");
		printf("m = %4d, ",         m);
		printf("n = %4d, ",         n);
		printf("k = %4d, ",         k);
		printf("nb = %4d, ",       nb);
		printf(" \n");
		printf("time = %16.8f (seconds)    performance = %10.3f (GFlop/sec) ", elapsed, perform);	
		printf(" \n");

	} else {

		printf("%6d %6d %6d %6d %16.8f %10.3f ", m, n, k, nb, elapsed, perform);

	} 

	if ( testing ){

		check_qq_orth( &orth, m, n, Q, m );
		if ( verbose ) printf("qq_orth  = %5.1e  \n",orth); else printf(" %5.1e  ",orth); 

		check_qr_repres( &repres, m, k, A, lda, Q, ldq, R, ldr );
		if ( verbose ) printf("qr_repres = %5.1e  \n",repres); else printf(" %5.1e  ",repres); 

	}

	if ( !verbose ) printf("\n");		

	free( A );
	free( Q );
	free( R );
	free( T );

	return 0;

}
