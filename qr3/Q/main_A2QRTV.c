#include "qr3.h"

int main(int argc, char ** argv) {

	int i, lda, ldq, ldt, m, n, verbose, testing;
	double *A, *As, *Q, *T;
	double orth, repres;
	double elapsed_ref, perform_ref;
	struct timeval tp;
	
	srand(0);

    	m         = 27;
    	n         = 20;
	lda       = -1;
	ldq       = -1;
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
	if( ldq < 0 ) ldq = m;
	if( ldt < 0 ) ldt = n;

	A  = (double *) malloc( lda * n * sizeof(double));
	As = (double *) malloc( lda * n * sizeof(double));
	Q  = (double *) malloc( ldq * n * sizeof(double));
 	T  = (double *) malloc( ldt * n * sizeof(double));

 	for(i = 0; i < lda * n; i++)
		*(A + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

	for(i = 0; i < ldq * n; i++)
		*(Q + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

	for(i = 0; i < ldt * n; i++)
		*(T + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

	LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, n, A, lda, As, lda );

	gettimeofday(&tp, NULL);
	elapsed_ref=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

	qr3_dA2QRTV( m, n, A, lda, Q, ldq, T, ldt );

	gettimeofday(&tp, NULL);
	elapsed_ref+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

	perform_ref = ((double) flops_org2r( m, n, n )) / elapsed_ref / 1.0e+9 ;

	if ( verbose ){ 

		printf("ORGQR");
		printf("m = %4d, ",         m);
		printf("n = %4d, ",         n);
		printf(" \n");
		printf(" time = %f    GFlop/sec = %f ", elapsed_ref, perform_ref);	
		printf(" \n ");

	} else {

		printf("%6d %6d %16.8f %10.3f ", m, n, elapsed_ref, perform_ref);

	} 

	if ( testing ){

		qr3_test_qq_orth_1( &orth, m, n, Q, ldq );		
		if ( verbose ) printf("qq_orth  = %5.1e  \n ",orth); else printf(" %5.1e  ",orth); 

		qr3_test_qr_repres_1( &repres, m, n, As, lda, Q, lda, A, lda );
		if ( verbose ) printf("qr_repres = %5.1e  \n ",repres); else printf(" %5.1e  ",repres); 

	}

	if ( !verbose ) printf("\n");		

	free( A  );
	free( As );
	free( Q  );
	free( T  );

	return 0;

}
