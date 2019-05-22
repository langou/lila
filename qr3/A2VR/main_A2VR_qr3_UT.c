#include "A2VR.h"

int main(int argc, char ** argv) {

	int i, lda, ldr, m, n, lwork, verbose, testing;
	double *A, *As, *R, *tau, *work;
	double orth, repres;
	double elapsed_ref, perform_ref;
	struct timeval tp;
	
	srand(0);

    	m         = 27;
    	n         = 20;
	lda       = -1;
	ldr       = -1;
	verbose   = 0;
	testing   = 1;

	for(i = 1; i < argc; i++){
		if( strcmp( *(argv + i), "-lda") == 0) {
			lda = atoi( *(argv + i + 1) );
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
	}

	if (( m < n )) { printf("\n\n We need n <= m\n\n"); return 0; }

	if( lda < 0 ) lda = m;
	if( ldr < 0 ) ldr = n;

	A  = (double *) malloc( lda * n * sizeof(double));
	As = (double *) malloc( lda * n * sizeof(double));
 	R  = (double *) malloc( ldr * n * sizeof(double));

 	for(i = 0; i < lda * n; i++)
		*(A + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

	for(i = 0; i < ldr * n; i++)
		*(R + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

	tau = (double *) malloc( n * sizeof(double));

	work = (double *) malloc( 1 * sizeof(double));
	LAPACKE_dgeqrf_work( LAPACK_COL_MAJOR, m, n, A, lda, tau, work, -1 ); 
	lwork = ((int) work[0]);
	LAPACKE_dorgqr_work( LAPACK_COL_MAJOR, m, n, n, A, lda, tau, work, -1 );
	if (lwork < ((int) work[0])) lwork = ((int) work[0]); 
	free( work );
	work = (double *) malloc( lwork * sizeof(double));

	LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, n, A, lda, As, lda );

	gettimeofday(&tp, NULL);
	elapsed_ref=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

	dgeqr3R_UT( m, n, A, lda, A, lda, R, ldr );

	gettimeofday(&tp, NULL);
	elapsed_ref+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

	long int flops_geqr3, int_m, int_n;
	int_m = m; int_n = n;
	flops_geqr3 = ( (( long int ) 18 ) * int_m * int_n * int_n - (( long int ) 5 ) * int_n * int_n * int_n ) / (( long int ) 6 );
	perform_ref = ( ((double) flops_geqr3 ) ) / elapsed_ref / 1.0e+9 ;

	if ( verbose ){ 

		printf("UT GEQR3");
		printf("m = %4d, ", m);
		printf("n = %4d, ", n);
		printf(" \n");
		printf(" time = %f    GFlop/sec = %f ", elapsed_ref, perform_ref);	
		printf(" \n ");

	} else {

		printf("%6d %6d %16.8f %10.3f ", m, n, elapsed_ref, perform_ref);

	} 

	if ( testing ){

		double *Q, *T;
		Q  = (double *) malloc( m * n * sizeof(double));
		T  = (double *) malloc( n * n * sizeof(double));
		LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'L', m, n, A, lda, Q, m );
		//LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'U', n, n, A, lda, T, n );
		//for(i=0;i<n;i++) tau[i] = T[i+i*n];
		dV2tau( m, n, A, lda, tau );
		LAPACKE_dlarft_work( LAPACK_COL_MAJOR, 'F', 'C', m, n, A, lda, tau, T, n );

		qr3_dorgqr( m, n, Q, m, T, n, tau );

		A2VR_test_qq_orth_1( &orth, m, n, Q, m );
		if ( verbose ) printf("qq_orth  = %5.1e  \n ",orth); else printf(" %5.1e  ",orth); 

		A2VR_test_qr_repres_1( &repres, m, n, As, lda, Q, m, R, ldr );
		if ( verbose ) printf("qr_repres = %5.1e  \n ",repres); else printf(" %5.1e  ",repres); 

		free( Q );
		free( T );

	}

	if ( !verbose ) printf("\n");		

	free( A  );
	free( As );
	free( R  );
	free( work );
	free( tau );

	return 0;

}
