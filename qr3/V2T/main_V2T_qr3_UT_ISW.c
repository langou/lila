#include "../src/qr2.h"
#include "../check/check.h"
#include "../flops/flops.h"

int main(int argc, char ** argv) {

	int i, lda, ldt, lwork, m, n, verbose, testing;
	double *A, *As, *T, *tau, *work;
	double orth, repres;
	double elapsed_ref, perform_ref;
	struct timeval tp;
	
	srand(0);

    	m         = 27;
    	n         = 20;
	lda       = -1;
	ldt       = -1;
	verbose   = 0;
	testing   = 1;

	for(i = 1; i < argc; i++){
		if( strcmp( *(argv + i), "-lda") == 0) {
			lda = atoi( *(argv + i + 1) );
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
	}

	if (( m < n )) { printf("\n\n We need n <= m\n\n"); return 0; }

	if( lda < 0 ) lda = m;
	if( ldt < 0 ) ldt = n;

	A  = (double *) malloc( lda * n * sizeof(double));
	As = (double *) malloc( lda * n * sizeof(double));
 	T  = (double *) malloc( ldt * n * sizeof(double));

 	for(i = 0; i < lda * n; i++)
		*(A + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

	for(i = 0; i < ldt * n; i++)
		*(T + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

	LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, n, A, lda, As, lda );

 	tau  = (double *) malloc( n * sizeof(double));
	work = (double *) malloc( 1 * sizeof(double));
	LAPACKE_dgeqrf_work( LAPACK_COL_MAJOR, m, n, A, lda, tau, work, -1 ); 
	lwork = ((int) work[0]);
	LAPACKE_dorgqr_work( LAPACK_COL_MAJOR, m, n, n, A, lda, tau, work, -1 );
	if (lwork < ((int) work[0])) lwork = ((int) work[0]); 
	free( work );
	work = (double *) malloc( lwork * sizeof(double));

	LAPACKE_dgeqrf_work( LAPACK_COL_MAJOR, m, n, A, lda, tau, work, lwork ); 

	gettimeofday(&tp, NULL);
	elapsed_ref=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

	qr2_larft_UT_ISW( m, n, A, lda, T, ldt, tau );

	gettimeofday(&tp, NULL);
	elapsed_ref+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

	long int flops;
	flops = flops_larft( m, n );
	perform_ref = ( ((double) flops ) ) / elapsed_ref / 1.0e+9 ;

	if ( verbose ){ 

		printf("QR3 UT ISW LARFT");
		printf("m = %4d, ", m);
		printf("n = %4d, ", n);
		printf(" \n");
		printf(" time = %f    GFlop/sec = %f ", elapsed_ref, perform_ref);	
		printf(" \n ");

	} else {

		printf("%6d %6d %16.8f %10.3f ", m, n, elapsed_ref, perform_ref);

	} 

	if ( testing ){

		double *Q;
		Q  = (double *) malloc( m * n * sizeof(double));
		LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'L', m, n, A, lda, Q, m );
		qr2_dorgqr_UT( m, n, Q, m, T, ldt, tau );

		check_qq_orth( &orth, m, n, Q, m );
		if ( verbose ) printf("qq_orth  = %5.1e  \n ",orth); else printf(" %5.1e  ",orth); 

		check_qr_repres( &repres, m, n, As, lda, Q, m, A, lda );
		if ( verbose ) printf("qr_repres = %5.1e  \n ",repres); else printf(" %5.1e  ",repres); 

		free( Q );
	}

	if ( !verbose ) printf("\n");		

	free( A  );
	free( As );
	free( T  );
	free( tau );
	free( work );

	return 0;

}
