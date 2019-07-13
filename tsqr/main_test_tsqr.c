#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include "cblas.h"
#include "lapacke.h"

extern int check_qq_orth( double *norm_orth_1, int m, int n, double *Q, int ldq  );
extern int check_qr_repres( double *norm_repres, int m, int n, double *A, int lda, double *Q, int ldq, double *R, int ldr );

extern int dlatsqr_( int *m, int *n, int *mb, int *nb, double *A, int *lda, double *T, int *ldt, double *work, int *lwork, int *info);

extern int dlamtsqr_( char *side, char *trans, int *m, int *n, int *k, int *mb, int *nb,
		double *A, int *lda, double *T, int *ldt, double *C, int *ldc, double *work, int *lwork, int *info );

int main(int argc, char ** argv) {

	int i, lda, ldq, ldr, ldt, m, n, verbose, testing;
	int lwork;
	double *A, *Q, *R, *T, *tau, *work, *Asave;
	double orth, repres;
	double elapsed, perform;
	struct timeval tp;
	
	srand(0);

    	m         = 27;
    	n         = 20;
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

	if ( m < n ) { printf("\n\n We need n <= m\n\n"); return 0; }

	if( lda < 0 ) lda = m;
	if( ldq < 0 ) ldq = m;
	if( ldr < 0 ) ldr = n;
	if( ldt < 0 ) ldt = n;

	A = (double *) malloc( lda * n * sizeof(double));
	Asave = (double *) malloc( lda * n * sizeof(double));
	Q = (double *) malloc( ldq * n * sizeof(double));
 	R = (double *) malloc( ldr * n * sizeof(double));

 	for(i = 0; i < lda * n; i++)
		*(A + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

 	for(i = 0; i < lda * n; i++)
		*(Asave + i) = *(A + i);

	for(i = 0; i < ldq * n; i++)
		*(Q + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

	for(i = 0; i < ldr * n; i++)
		*(R + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

//	tau = (double *) malloc( n * sizeof(double));

//	work = (double *) malloc( 1 * sizeof(double));
//	LAPACKE_dgeqrf_work( LAPACK_COL_MAJOR, m, n, Q, ldq, tau, work, -1 ); 
//	lwork = ((int) work[0]);
//	LAPACKE_dorgqr_work( LAPACK_COL_MAJOR, m, n, n, Q, ldq, tau, work, -1 );
//	if (lwork < ((int) work[0])) lwork = ((int) work[0]); 
//	free( work );
//	work = (double *) malloc( lwork * sizeof(double));

//	gettimeofday(&tp, NULL);
//	elapsed=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

//	LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, n, A, lda, Q, ldq );
//	LAPACKE_dgeqrf_work( LAPACK_COL_MAJOR, m, n, Q, ldq, tau, work, lwork ); 
//	LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'U', n, n, Q, ldq, R, ldr );
//	LAPACKE_dorgqr_work( LAPACK_COL_MAJOR, m, n, n, Q, ldq, tau, work, lwork );

//	gettimeofday(&tp, NULL);
//	elapsed+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

//	free( tau );
//	free( work );

	int mb = m;
	int nb = n / 2;
	int info;

 	T = (double *) malloc( ldt * n * ceil((m-n)/(mb-n)) * sizeof(double));

	lwork = -1;
	work = (double *) malloc( 1 * sizeof(double));
	dlatsqr_( &m, &n, &mb, &nb, A, &lda, T, &ldt, work, &lwork, &info );
	lwork = ((int) work[0]);
	free( work );
	work = (double *) malloc( lwork * sizeof(double));

	dlatsqr_( &m, &n, &mb, &nb, A, &lda, T, &ldt, work, &lwork, &info );

	LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'U', n, n, A, lda, R, ldr );

	LAPACKE_dlaset_work( LAPACK_COL_MAJOR, 'A', m, n, 0.0e+00, 1.0e+00, Q, ldq );

	dlamtsqr_( "L", "N", &m, &n, &n, &mb, &nb, A, &lda, T, &ldt, Q, &ldq, work, &lwork, &info );

	free( work );
	free( T );

//	perform = ((double) flops_lapack_org2r( m, n, k ) + (double) flops_lapack_geqr2( m, k ) ) / elapsed / 1.0e+9 ;
	perform = 0.0e+00 ;

	if ( verbose ){ 

		printf("ORGQR  ");
		printf("m = %4d, ",         m);
		printf("n = %4d, ",         n);
		printf(" \n");
		printf("time = %16.8f (seconds)    performance = %10.3f (GFlop/sec) ", elapsed, perform);	
		printf(" \n");

	} else {

		printf("%6d %6d %16.8f %10.3f ", m, n, elapsed, perform);

	} 

	if ( testing ){

		check_qq_orth( &orth, m, n, Q, ldq );
		if ( verbose ) printf("qq_orth  = %5.1e  \n",orth); else printf(" %5.1e  ",orth); 

		check_qr_repres( &repres, m, n, Asave, lda, Q, ldq, R, ldr );
		if ( verbose ) printf("qr_repres = %5.1e  \n",repres); else printf(" %5.1e  ",repres); 

	}

	if ( !verbose ) printf("\n");		

	free( Asave );
	free( A );
	free( Q );
	free( R );

	return 0;

}

int check_qr_repres( double *norm_repres, int m, int n, double *A, int lda, double *Q, int ldq, double *R, int ldr ){

	double normA, *work;
	int ii, jj;

	normA = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', m, n, A, lda, NULL );

	work  = (double *) malloc(m * n * sizeof(double));
	LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, n, Q, ldq, work, m );
	cblas_dtrmm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, m, n, (1.0e+00), R, ldr, work, m );
 	for(ii = 0; ii < m; ii++) for(jj = 0; jj < n; jj++) work[ ii+jj*m ] -= A[ ii+jj*lda ];
	(*norm_repres) = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', m, n, work, m, NULL );
	free( work );

	(*norm_repres) = (*norm_repres) / normA;

	return 0;
}

int check_qq_orth( double *norm_orth_1, int m, int n, double *Q, int ldq  ){

	double *work;

	work  = (double *) malloc( n * n * sizeof(double));
	LAPACKE_dlaset( LAPACK_COL_MAJOR, 'A', n, n, (0e+00), (1e+00), work, n );
	cblas_dsyrk( CblasColMajor, CblasUpper, CblasTrans, n, m, 1.0e+00, Q, ldq, -1.0e+00, work, n );
	(*norm_orth_1) = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', n, n, work, n, NULL );
	free( work );

	return 0;
}
