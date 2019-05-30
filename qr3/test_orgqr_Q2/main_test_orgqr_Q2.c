#include "../src/qr2.h"
#include "../check/check.h"
#include "../flops/flops.h"
#include "null.h"

int main(int argc, char ** argv) {

	int i, lda, ldq1, ldq2, ldr, ldt, lwork, m, n, k, nb, verbose, testing;
	double *A, *Q1, *Q2, *R, *T, *tau, *work;
	double orth, repres;
	double elapsed, perform_ref;
	struct timeval tp;
	
	srand(0);

    	m         = 27;
    	n         = 20;
    	k         = 10;
    	nb        = 20;
	lda       = -1;
	ldq1      = -1;
	ldq2      = -1;
	ldr       = -1;
	ldt       = -1;
	verbose   =  0;
	testing   =  1;

	for(i = 1; i < argc; i++){
		if( strcmp( *(argv + i), "-lda") == 0) {
			lda = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-ldq1") == 0) {
			ldq1 = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-ldq2") == 0) {
			ldq2 = atoi( *(argv + i + 1) );
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

	if( lda  < 0 ) lda  = m;
	if( ldq1 < 0 ) ldq1 = m;
	if( ldq2 < 0 ) ldq2 = m;
	if( ldr  < 0 ) ldr  = k;
	if( ldt  < 0 ) ldt  = k;

	A  = (double *) malloc( lda  * n   * sizeof(double));
	Q1 = (double *) malloc( ldq1 * k   * sizeof(double));
	Q2 = (double *) malloc( ldq2 * n-k * sizeof(double));
 	R  = (double *) malloc( ldr  * k   * sizeof(double));
 	T  = (double *) malloc( ldt  * k   * sizeof(double));

 	for(i = 0; i < lda * n; i++)
		*(A + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

	for(i = 0; i < ldr  * k; i++)
		*(R + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

	lwork = n;
	work  = (double *) malloc( lwork * n * sizeof(double));
	tau   = (double *) malloc( k * sizeof(double));

	LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m,   k,       A, lda, Q1, ldq1 );
	LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, n-k, A+k*lda, lda, Q2, ldq2 );

	qr2_dgeqr3_R( m, k, Q1, ldq1, T, ldt, R, ldr );
	for(i=0;i<k;i++) tau[i] = T[i+i*ldt];

	gettimeofday(&tp, NULL);
	elapsed=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);


	double *Q;
	Q = (double *) malloc( m * n * sizeof(double));	
	LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, k, Q1, ldq1, Q, m );

//	The reason for having the variants,
//	I can not seem to write in Q2. It ruins the test checks if one of the following is uncommented 
//	except fpr when I call Q only, or call the /src/ file

//	lapack_our_dorgq2r( m, n, k, nb, Q1, ldq1, Q2, ldq2, tau, work, lwork );
//	lapack_our_dorgq2r( m, n, k, nb, Q, m, Q2, ldq2, tau, work, lwork );
//	lapack_our_dorgq2r( m, n, k, nb, Q, m, Q+k*m, m, tau, work, lwork );

//	qr3_aux_dorgq2r( m, n, k, Q1, ldq1, T, ldt, Q2, ldq2 );
//	qr3_aux_dorgq2r( m, n, k, Q, m, T, ldt, Q2, ldq2 );
//	qr3_aux_dorgq2r( m, n, k, Q, m, T, ldt, Q+k*m, m );

//	/src -- files
//	lapack_our_dorgqr_Q2( m, n, k, nb, Q, m, tau, work, lwork );



	gettimeofday(&tp, NULL);
	elapsed+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

	qr2_dorgqr3( m, k, Q1, ldq1, T, ldt, tau );

	free( tau );
	free( work );

	perform_ref = ( ((double) flops_lapack_org2r( m, n, k ) ) + ((double) flops_lapack_geqr2( m, n ) ) ) / elapsed / 1.0e+9 ;

	if ( verbose ){ 

		printf("ORGQR and QRGQ2R");
		printf("m = %4d, ",           m);
		printf("n = %4d, ",           n);
		printf("k = %4d, ",           k);
		printf("nb = %4d, ",         nb);
		printf(" \n");
		printf(" time = %f    GFlop/sec = %f ", elapsed, perform_ref);	
		printf(" \n ");

	} else {

		printf("%6d %6d %6d %6d %16.8f %10.3f ", m, n, k, nb, elapsed, perform_ref);

	} 

	if ( testing ){

		check_qq_orth( &orth, m, k, Q1, ldq1 );		
		//check_qq_orth( &orth, m, k, Q, m );		
		if ( verbose ) printf("q1q1_orth  = %5.1e  \n ",orth); else printf(" %5.1e  ",orth);

		check_qr_repres( &repres, m, k, A, lda, Q1, ldq1, R, ldr );
		//check_qr_repres( &repres, m, k, A, lda, Q, m, R, ldr );
		if ( verbose ) printf("qr_repres = %5.1e  \n ",repres); else printf(" %5.1e  ",repres); 


	//	double normA, norm_repres_q2;
	//	normA = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', m, n, A, lda, NULL );
	//	work  = (double *) malloc(m * n * sizeof(double));
	//	LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m,   k, Q1, ldq1, work, m );
	//	LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, n-k, Q2, ldq2, work+k*m, m );
	//	cblas_dtrmm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, m, k, (+1.0e+00), R, ldr, work, m );
 	//	{ int ii,jj; for(ii = 0; ii < m; ii++) for(jj = 0; jj < k; jj++) work[ ii+jj*m ] -= A[ ii+jj*lda ]; }
	//	cblas_dgemm( CblasColMajor, CblasTrans, CblasNoTrans, n-k, n-k, m, (+1.0e+00), work+k*m, m, A+k*lda, lda, (+0.0e+00), work+k*m, m );
	//	//cblas_dgemm( CblasColMajor, CblasTrans, CblasNoTrans, n-k, n, m, (+1.0e+00), work+k*m, m, A, lda, (+0.0e+00), work+k*m, m );
	//	LAPACKE_dlaset( LAPACK_COL_MAJOR, 'A', m-n, n, (+0.0e00), (+0.0e00), work+(m-n+k)+(n-k)*m, m );
	//	norm_repres_q2 = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', m, n, work, m, NULL );
	//	free( work );
	//	norm_repres_q2 = norm_repres_q2 / normA;
	//	if ( verbose ) printf("q_repres  = %5.1e  \n ",norm_repres_q2); else printf(" %5.1e  ",norm_repres_q2); 


	}

	if ( !verbose ) printf("\n");		

	free( A  );
	free( Q1 );
	free( Q2 );
	free( R  );
	free( Q  );
	free( T  );

	return 0;

}
