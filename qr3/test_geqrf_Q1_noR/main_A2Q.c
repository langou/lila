#include "../src/qr2.h"
#include "../check/check.h"
#include "../flops/flops.h"

int main(int argc, char ** argv) {

	int i, lda, ldq, m, n, verbose, testing;
	double *A, *Q, *tau;
	double orth, repres;
	double elapsed, perform_rel, perform_abs;
	struct timeval tp;
	
	srand(0);

    	m         = 27;
    	n         = 20;
	lda       = -1;
	ldq       = -1;
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
	}

	if (( m < n )) { printf("\n\n We need n <= m\n\n"); return 0; }

	if( lda < 0 ) lda = m;
	if( ldq < 0 ) ldq = m;

	A  = (double *) malloc( lda * n * sizeof(double));
	Q  = (double *) malloc( ldq * n * sizeof(double));

 	for(i = 0; i < lda * n; i++)
		*(A + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

 	for(i = 0; i < ldq * n; i++)
		*(Q + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

	tau = (double *) malloc( n * sizeof(double));

	LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, n, A, lda,  Q, ldq );

	gettimeofday(&tp, NULL);
	elapsed=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

	dgeqr3_Q( m, n, Q, ldq, Q, ldq );
	for(i=0;i<n;i++) tau[i] = Q[i+i*ldq];
	qr2_dorgqr3( m, n, Q, ldq, Q, ldq, tau );

	gettimeofday(&tp, NULL);
	elapsed+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

	perform_rel = ( ((double) flops_lapack_geqr2( m, n ) ) + ((double) flops_lapack_org2r( m, n, n ) ) ) / elapsed / 1.0e+9 ;

	perform_abs = ( ((double) flops_geqr3( m, n ) ) + ((double) flops_lapack_org2r( m, n, n ) ) ) / elapsed / 1.0e+9 ;

	if ( verbose ){ 

		printf("GEQRF ORGQR NO R");
		printf("m = %4d, ", m);
		printf("n = %4d, ", n);
		printf(" \n");
		printf(" time = %f    GFlop/sec (rel) = %f GFlop/sec (abs) = %f ", elapsed, perform_rel, perform_abs);	
		printf(" \n ");

	} else {

		printf("%6d %6d %16.8f %10.3f %10.3f ", m, n, elapsed, perform_rel, perform_abs);

	} 

	if ( testing ){

		double *R, *tmp;
		R  = (double *) malloc( n * n * sizeof(double));
		tmp  = (double *) malloc( m * n * sizeof(double));
		LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, n, A, lda, tmp, m );
		qr2_dgeqr3_R( m, n, tmp, m, tmp, m, R, n );

		check_qq_orth( &orth, m, n, Q, ldq );
		if ( verbose ) printf("qq_orth  = %5.1e  \n ",orth); else printf(" %5.1e  ",orth); 

		check_qr_repres( &repres, m, n, A, lda, Q, ldq, R, n );
		if ( verbose ) printf("qr_repres = %5.1e  \n ",repres); else printf(" %5.1e  ",repres); 

		free( tmp );
		free( R );

	}

	if ( !verbose ) printf("\n");		

	free( A  );
	free( Q  );
	free( tau );

	return 0;

}
