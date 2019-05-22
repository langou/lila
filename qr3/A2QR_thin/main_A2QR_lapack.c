#include "../src/qr2.h"
#include "../check/check.h"

int main(int argc, char ** argv) {

	int i, lda, ldq, ldr, m, n, lwork, verbose, testing;
	double *A, *Q, *R, *tau, *work;
	double orth, repres;
	double elapsed_ref, perform_ref;
	struct timeval tp;
	
	srand(0);

    	m         = 27;
    	n         = 20;
	lda       = -1;
	ldq       = -1;
	ldr       = -1;
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
	}

	if (( m < n )) { printf("\n\n We need n <= m\n\n"); return 0; }

	if( lda < 0 ) lda = m;
	if( ldq < 0 ) ldq = m;
	if( ldr < 0 ) ldr = n;

	A  = (double *) malloc( lda * n * sizeof(double));
	Q  = (double *) malloc( ldq * n * sizeof(double));
 	R  = (double *) malloc( ldr * n * sizeof(double));

 	for(i = 0; i < lda * n; i++)
		*(A + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

 	for(i = 0; i < ldq * n; i++)
		*(Q + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

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

	LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, n, A, lda,  Q, ldq );

	gettimeofday(&tp, NULL);
	elapsed_ref=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

	LAPACKE_dgeqrf_work( LAPACK_COL_MAJOR, m, n, Q, ldq, tau, work, lwork ); 
	LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'U', n, n, Q, ldq, R, ldr );
	LAPACKE_dorgqr_work( LAPACK_COL_MAJOR, m, n, n, Q, ldq, tau, work, lwork );

	gettimeofday(&tp, NULL);
	elapsed_ref+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

	long int flops_geqrf, flops_orgqr, int_m, int_n;
	int_m = m; int_n = n;
	flops_geqrf =  ( (( long int ) 6 ) * int_m * int_n * int_n - (( long int ) 2 ) * int_n * int_n * int_n + (( long int ) 3 ) * int_m * int_n 
	+ (( long int ) 3 ) * int_n * int_n + (( long int ) 42 ) * int_n ) / (( long int ) 3 );
	flops_orgqr = ( (( long int ) 12 ) * int_m * int_n * int_n - (( long int ) 6 ) * ( int_m + int_n ) * int_n * int_n + (( long int ) 4 ) * int_n * int_n * int_n 
	+ (( long int ) 9 ) * int_n * int_n - (( long int ) 3 ) * int_m * int_n - (( long int ) 3 ) * int_n * int_n - (( long int ) 4 ) * int_n ) / (( long int ) 3 );
	perform_ref = ( ((double) flops_geqrf ) + ((double) flops_orgqr ) ) / elapsed_ref / 1.0e+9 ;

	if ( verbose ){ 

		printf("LAPACK GEQRF ORGQR");
		printf("m = %4d, ", m);
		printf("n = %4d, ", n);
		printf(" \n");
		printf(" time = %f    GFlop/sec = %f ", elapsed_ref, perform_ref);	
		printf(" \n ");

	} else {

		printf("%6d %6d %16.8f %10.3f ", m, n, elapsed_ref, perform_ref);

	} 

	if ( testing ){

		check_qq_orth( &orth, m, n, Q, ldq );
		if ( verbose ) printf("qq_orth  = %5.1e  \n ",orth); else printf(" %5.1e  ",orth); 

		check_qr_repres( &repres, m, n, A, lda, Q, ldq, R, ldr );
		if ( verbose ) printf("qr_repres = %5.1e  \n ",repres); else printf(" %5.1e  ",repres); 

	}

	if ( !verbose ) printf("\n");		

	free( A  );
	free( Q  );
	free( R  );
	free( work );
	free( tau );

	return 0;

}
