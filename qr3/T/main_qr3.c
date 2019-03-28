#include "lila.h"

int main(int argc, char ** argv) {

	int i, info, lda, ldt, m, n, nx, verbose, testing;
	int lwork, leaf, vrtq, *lila_param;
	double *A, *As, *T, *work=NULL;
	double elapsed_ref1, perform_ref1, elapsed_ref2, perform_ref2;
	struct timeval tp;

	srand(0);

    	m         = 37;
    	n         = 23;
	lda       = -1;
	nx        = 100;
	leaf      = 0;
	verbose   = 0;
	testing   = 1;
	vrtq      = 3;


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
		if( strcmp( *(argv + i), "-nx") == 0) {
			nx  = atoi( *(argv + i + 1) );
			i++;
		}
	}

	if( m < n ){ printf("\n\n YOUR CHOICE OF n AND ii HAVE MADE YOU LARGER THAN m, PLEASE RECONSIDER \n\n"); return 0; }

	if( lda < 0 ) lda = m;
	    	      ldt = n;

	lila_param = (int *) malloc(6 * sizeof(int));
	lila_param[ 0 ] = 0;
	lila_param[ 1 ] = leaf;
	lila_param[ 2 ] = -1;
	lila_param[ 3 ] = nx;
	lila_param[ 4 ] = vrtq;

	A  = (double *) malloc( lda * (n) * sizeof(double));
 	As = (double *) malloc( lda * (n) * sizeof(double));
	T  = (double *) malloc( ldt * (n) * sizeof(double));

 	for(i = 0; i < lda * (n); i++)
		*(A + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, n, A, lda, As, lda );

	lwork = lila_query_dgeqrf_w03_recursive( lila_param, m, n, 0, n, A, lda, T, ldt, NULL, -1, work=NULL, -1 );
	work  = (double *) malloc( lwork * sizeof(double));

	gettimeofday(&tp, NULL);
	elapsed_ref1=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

	// Getting V for the checks below and the consrtuction of T
	lila_param[4] = 1; //   ==>  VRT
	lila_dgeqrf_recursive( lila_param, m, n, 0, n, A, lda, T, ldt, NULL, -1, work, lwork ); 
	gettimeofday(&tp, NULL);
	elapsed_ref1+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
	perform_ref1 = ( 2.0e+00 * ((double) m) * ((double) n) * ((double) n) - 2.0e+00 / 3.0e+00 * ((double) n) * ((double) n) * ((double) n) )  / elapsed_ref1 / 1.0e+9 ;
	info = LAPACKE_dlaset( LAPACK_COL_MAJOR, 'A', n, n, (0e+00), (0e+00), T, ldt );


	gettimeofday(&tp, NULL);
	elapsed_ref2=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

	lila_param[4] = 3; //   ==> T
	lila_dgeqrf_recursive( lila_param, m, n, 0, n, A, lda, T, ldt, NULL, -1, work, lwork ); 

	gettimeofday(&tp, NULL);
	elapsed_ref2+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
	perform_ref2 = ( (2/3)*((double) n)*((double) n)*((double) n) )  / elapsed_ref2 / 1.0e+9 ;

	free( work );

	if ( verbose ){ 

		printf("larft_recursive - ");
		printf("m = %4d, ",         m);
		printf("n = %4d, ",         n);
		printf("lda = %4d, ",     lda);
		printf("nx = %4d ", nx); 
		printf(" \n");
		printf(" time = %f    GFlop/sec = %f ", elapsed_ref1, perform_ref1);	
		printf(" timeT = %f    GFlop/secT = %f ", elapsed_ref2, perform_ref2);	
		printf(" \n ");

	} else {

		printf("%6d %6d %6d %6d %16.8f %10.3f %16.8f %10.3f ", m, n, lda, nx, elapsed_ref1, perform_ref1, elapsed_ref2, perform_ref2);

	} 

	if ( testing ){

		double repres_1, repres_2, repres_3;
	
//		vrtq == 0 || 2    ===>   we have Q to check
//		vrtq == 1 || 3    ===>   we don't
	
//		|| (apply H) A - [R;0] || / || A ||
		repres_1 = lila_test_r_repres_2( lila_param, m, n, 0, n, As, lda, T, ldt, A, lda );

//		create m-by-m H, then || I - HH^T ||; || H A - [R;0] || / || A ||
		repres_2 = lila_test_hh_repres( lila_param, m, n, 0, n, As, lda, T, ldt, A, lda );

//		create m-by-n Q, then || I - Q^T Q ||; || A - Q R || / || A ||
		repres_3 = lila_test_vt_repres( lila_param, m, n, 0, n, As, lda, T, ldt, A, lda );

		if ( verbose ) printf("r_repres  = %5.1e  \n ",repres_1); else printf(" %5.1e  ",repres_1); 
		if ( verbose ) printf("h_q_orth  = %5.1e  \n ",repres_2); else printf(" %5.1e  ",repres_2); 
		if ( verbose ) printf("vt_repres = %5.1e  \n ",repres_3); else printf(" %5.1e  ",repres_3); 

	}

	if ( !verbose ) printf("\n");		

	free( A );
	free( As );
	free( T );
	free( lila_param );

	return 0;
}
