#include "lila.h"

int main(int argc, char ** argv) {

	int i, info, lda, ldt, m, n, verbose, testing, lwork;
	double *A, *As, *T, *work=NULL;
	double elapsed_ref1, perform_ref1;
	struct timeval tp;

	srand(0);

    	m         = 87;
    	n         = 23;
	lda       = -1;
	verbose   = 0;
	testing   = 1;


	for(i = 1; i < argc; i++){
		if( strcmp( *(argv + i), "-lda") == 0) {
			lda  = atoi( *(argv + i + 1) );
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

	if( m < n ){ printf("\n\n YOUR CHOICE OF n AND ii HAVE MADE YOU LARGER THAN m, PLEASE RECONSIDER \n\n"); return 0; }

	if( lda < 0 ) lda = m;
	    	      ldt = n;

        A  = (double *) malloc(lda * (n) * sizeof(double));
	As = (double *) malloc(lda * (n) * sizeof(double));
   	T  = (double *) malloc(ldt * (n) * sizeof(double));

 	for(i = 0; i < lda * (n); i++)
		*(A + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, n, A, lda, As, lda );

	lwork = n*n;
	work  = (double *) malloc( lwork * sizeof(double));

	gettimeofday(&tp, NULL);
	elapsed_ref1=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

	info = dgeqr3( m, n, A, lda, T, ldt ); 

	gettimeofday(&tp, NULL);
	elapsed_ref1+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
	perform_ref1 = ( 4.0e+00 * ((double) m) * ((double) n) * ((double) n) - 4.0e+00 / 3.0e+00 * ((double) n) * ((double) n) * ((double) n) )  / elapsed_ref1 / 1.0e+9 ;

	free( work );


	if ( verbose ){ 

		printf("dgeqrf_LAPACK - ");
		printf("m = %4d, ",         m);
		printf("n = %4d, ",         n);
		printf("lda = %4d, ",     lda);
		printf(" \n");
		printf(" time = %f    GFlop/sec = %f ", elapsed_ref1, perform_ref1);	
		printf(" \n ");

	} else {

		printf("%6d %6d %6d %16.8f %10.3f ", m, n, lda, elapsed_ref1, perform_ref1);

	} 

	if ( testing ){

		double repres_1, repres_2, repres_3;
		int *lila_param;
	        lila_param  = (int *) malloc(5 * sizeof(int));
		for( i = 0; i < 6; i++) lila_param[i] = 1;
	
//		|| (apply H) A - [R;0] || / || A ||
		repres_1 = lila_test_r_repres_2( lila_param, m, n, 0, n, As, lda, T, ldt, A, lda );

//		create m-by-m H, then || I - HH^T ||; || H A - [R;0] || / || A ||
		repres_2 = lila_test_hh_repres( lila_param, m, n, 0, n, As, lda, T, ldt, A, lda );

//		create m-by-n Q, then || I - Q^T Q ||; || A - Q R || / || A ||
		repres_3 = lila_test_vt_repres( lila_param, m, n, 0, n, As, lda, T, ldt, A, lda );

		free( lila_param );

		if ( verbose ) printf("r_repres  = %5.1e  \n ",repres_1); else printf(" %5.1e  ",repres_1); 
		if ( verbose ) printf("h_q_orth  = %5.1e  \n ",repres_2); else printf(" %5.1e  ",repres_2); 
		if ( verbose ) printf("vt_repres = %5.1e  \n ",repres_3); else printf(" %5.1e  ",repres_3); 

	}

	if ( !verbose ) printf("\n");		

	free( A );
	free( As );
	free( T );

	return 0;
}
