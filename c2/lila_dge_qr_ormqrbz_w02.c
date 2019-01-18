#include "lila.h"

int lila_dge_qr_ormqrbz_w02( int m, int n, int k, int i, int j, int mt, double *A, int lda, double *Q, int ldq, double *T, int ldt, double *work, int lwork ){

	double *Aii, *Qij;
	int ml;
	double *Tii;
	int ldwork;
	
	Aii = A + i + i*lda;
	Qij = Q + i + j*ldq;
	ml = m - i;

	ldwork = k;

	Tii = T + i + i*ldt;

 	cblas_dgemm( CblasColMajor, CblasTrans, CblasNoTrans, k, n, m-k-i, (1.0e+00), Aii+k, lda, Qij+k, lda, (0.0e+00), work, ldwork );

	cblas_dtrmm( CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, k, n, (1.0e+00), Tii, ldt, work, ldwork );

	LAPACKE_dlacpy( LAPACK_COL_MAJOR, 'A', k, n, work, ldwork, Qij, ldq );

	cblas_dtrmm( CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, k, n, (-1.0e+00), Aii, lda, Qij, ldq );

	cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, m-k-i, n, k, (-1.0e+00), Aii+k, lda, work, ldwork, (1.0e+00), Qij+k, ldq );

	return 0;

}


/*
	int info;
	double *V;
	double normV_square;
	int vl;
	double *tau;
	tau = (double *) malloc( k * sizeof(double));

	for( jj = 0, vl=ml-1, V = Aii+1; jj < k; jj++, vl--, V+=(lda+1) ){
		normV_square = ( 1.0e+00 ) + cblas_ddot( vl, V, 1, V, 1 );
		tau[jj] = ( 2.0e+00 ) / normV_square ;
	}

	info = LAPACKE_dlaset( LAPACK_COL_MAJOR, 'A', k, n, (0e+00), (0e+00), Qij, ldq );
	info = LAPACKE_dormqr_work( LAPACK_COL_MAJOR, 'L', 'N', ml, n, k, Aii, lda, tau, Qij, ldq, work, lwork );
	free( tau );
*/

/*
	int info;
	double *V;
	double normV_square;
	int vl;
	double *tau;
	tau = (double *) malloc( k * sizeof(double));

	for( jj = 0, vl=ml-1, V = Aii+1; jj < k; jj++, vl--, V+=(lda+1) ){
		normV_square = ( 1.0e+00 ) + cblas_ddot( vl, V, 1, V, 1 );
		tau[jj] = ( 2.0e+00 ) / normV_square ;
	}

	info = LAPACKE_dlaset( LAPACK_COL_MAJOR, 'A', k, n, (0e+00), (0e+00), Qij, ldq );

	double *Tii;
	int ldwork;
	ldwork = k;
	Tii = T + i + i*ldt;
	info = LAPACKE_dlarft_work ( LAPACK_COL_MAJOR, 'F', 'C', ml, k, Aii, lda, tau, Tii, ldt);
	info = LAPACKE_dlarfb_work ( LAPACK_COL_MAJOR, 'L', 'N', 'F', 'C', ml, n, k, Aii, lda, Tii, ldt, Qij, ldq, work, ldwork );
	free( tau );
*/


