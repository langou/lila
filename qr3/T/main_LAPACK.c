#include "lila.h"

extern int xN2T( int n, double *tau, double *T, int ldt );
extern int xV2N( int n, double *T, int ldt );



int main(int argc, char ** argv) {

	int i, j, info, lda, ldt, m, n, nx, verbose, testing;
	int lwork, leaf, vrtq, *lila_param;
	double *A, *As, *T, *tau, *work=NULL;
	double elapsed_ref1, perform_ref1, elapsed_ref2, perform_ref2;
	struct timeval tp;

	srand(0);

    	m         = 27;
    	n         = 17;
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

	tau   = (double *) malloc( n * sizeof(double));
	lwork = -1;
	work  = (double *) malloc( 1 * sizeof(double));
	info = LAPACKE_dgeqrf_work( LAPACK_COL_MAJOR, m, n, A, lda, tau, work, lwork ); 
	lwork = (int) work[0];
	free(work);
	work  = (double *) malloc( lwork * sizeof(double));

	gettimeofday(&tp, NULL);
	elapsed_ref1=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
	info = LAPACKE_dgeqrf_work( LAPACK_COL_MAJOR, m, n, A, lda, tau, work, lwork ); 
	gettimeofday(&tp, NULL);
	elapsed_ref1+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
	perform_ref1 = ( 2.0e+00 * ((double) m) * ((double) n) * ((double) n) - 2.0e+00 / 3.0e+00 * ((double) n) * ((double) n) * ((double) n) )  / elapsed_ref1 / 1.0e+9 ;

	gettimeofday(&tp, NULL);
	elapsed_ref2=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

	//info = LAPACKE_dlarft_work( LAPACK_COL_MAJOR, 'F', 'C', m, n, A, lda, tau, T, ldt);

	//nV2T( m, n, A, lda, tau, T, ldt, work, lwork );

	for(i=0;i<n;i++){T[i+i*n] = 1.0e+00; for(j=0;j<i;j++){ T[j+i*n] = A[i+j*lda];}}
	xV2N( n, T, ldt );
	cblas_dsyrk( CblasColMajor, CblasUpper, CblasTrans, n, m-n, (+1.0e+00), A+n, lda, (+1.0e+00), T, ldt );
	if (m==n) T[n-1+ldt*(n-1)]=0.0e+00;
	xN2T( n, tau, T, ldt );

	gettimeofday(&tp, NULL);
	elapsed_ref2+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);
	perform_ref2 = ( ((double) m)*((double) n)*((double) n)-(1/3)*((double) n)*((double) n)*((double) n) )  / elapsed_ref2 / 1.0e+9 ;

	//    CHECKS 
//	for(i=0;i<n;i++){ printf(" %+6.4f ", tau[i] ); }printf("\n");
	//T11
//	for(i=0;i<n1;i++){for(j=0;j<n1;j++){ printf(" %+1.2f ", Tc[i+j*ldt]-T[i+j*ldt]); }printf("\n");}
//	printf("\n");
	//T22
//	for(i=n1;i<n1+n2;i++){for(j=n1;j<n1+n2;j++){ printf(" %+1.2f ", Tc[i+j*ldt]-T[i+j*ldt]); }printf("\n");}
//	printf("\n");
	//T12
//	for(i=0;i<n1;i++){for(j=n1;j<n1+n2;j++){ printf(" %+1.2f ", Tc[i+j*ldt]-T[i+j*ldt]); }printf("\n");}
//	printf("\n");

//	free( Tc );

	free( work );
	free( tau );

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

		double *Q, orth;
		double repres_1;
		int ldq;
		ldq = m;
		Q  = (double *) malloc( ldq * n * sizeof(double));
		work  = (double *) malloc( n * n * sizeof(double));

		info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'L', m, n, A, lda, Q, ldq );
		info = LAPACKE_dlaset( LAPACK_COL_MAJOR, 'U', n, n, (0.0e+00), (1.0e+00), Q, ldq);
		info = LAPACKE_dlaset( LAPACK_COL_MAJOR, 'A', n, n, (0.0e+00), (0.0e+00), work, n);
		info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'U', n, n, T, ldt, work, n );
//		note that this is a triangle times triangle so we could have dtrtrmm() and dived # of FLOPS by 3x
		cblas_dtrmm( CblasColMajor, CblasRight, CblasLower, CblasTrans, CblasUnit, n, n, (+1.0e+00), A, lda, work, n );
//		note that the top n-by-n part is a lower times a upper and can be done dlauum
		cblas_dtrmm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, m, n, (-1.0e+00), work, n, Q, ldq );
	 	for(i = 0; i < n; i++) Q[ i + ldq * i ] = 1.00e+00 + Q[ i + ldq * i ];
//		|| I - Q^T Q ||
		orth = lila_test_qq_orth_1( m, n, 0, Q, ldq );		
//		create m-by-n Q, then || I - Q^T Q ||; || A - Q R || / || A ||
		repres_1 = lila_test_vt_repres( lila_param, m, n, 0, n, As, lda, T, ldt, A, lda );
		if ( verbose ) printf("qq_orth  = %5.1e  \n ",orth); else printf(" %5.1e  ",orth); 
		if ( verbose ) printf("vt_repres = %5.1e  \n ",repres_1); else printf(" %5.1e  ",repres_1); 

		free( work );
		free( Q );

	}

	if ( !verbose ) printf("\n");		

	free( A );
	free( As );
	free( T );
	free( lila_param );

	return 0;
}
