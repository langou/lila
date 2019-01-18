#include "lila.h"

int lila_dge_qr_ormqrf_w00( int m, int n, int k, int i, int j, int mt, double *A, int lda, double *T, int ldt, double *work, int lwork ){

	int info, vl, ml, jj;
	double normV_square;
	double *V, *tau, *Aii, *Tii, *Aij;

	Aii = A + i + i*lda;
	Aij = A + i + j*lda;
	Tii = T + i + i*ldt;
	ml = m - i;

	tau = (double *) malloc( k * sizeof(double));

	for(jj = 0, vl=ml-1, V = Aii+1; jj < k; jj++,vl--,V+=(lda+1) ){
		normV_square = ( 1.0e+00 ) + cblas_ddot( vl, V, 1, V, 1 );
		tau[jj] = ( 2.0e+00 ) / normV_square ;
	}

	info = LAPACKE_dormqr_work( LAPACK_COL_MAJOR, 'L', 'T', ml, n, k, Aii, lda, tau, Aij, lda, work, lwork );

	free( tau );
	return 0;

}

/*
	int ldwork, info, ii, jj, vl, ml;
	double normV_square;
	double *V, *tau, *Aii, *Tii, *Aij;

	Aii = A + i + i*lda;
	Aij = A + i + j*lda;
	Tii = T + i + i*ldt;
	ml = m - i;
	ldwork = n;

	tau = (double *) malloc( k * sizeof(double));
	for(jj = 0, vl=ml-1, V = Aii+1; jj < k; jj++,vl--,V+=(lda+1) ){
		normV_square = ( 1.0e+00 ) + cblas_ddot( vl, V, 1, V, 1 );
		tau[jj] = ( 2.0e+00 ) / normV_square ;
	}

	info = LAPACKE_dlarft_work ( LAPACK_COL_MAJOR, 'F', 'C', ml, k, Aii, lda, tau, Tii, ldt);

	info = LAPACKE_dlarfb_work ( LAPACK_COL_MAJOR, 'L', 'T', 'F', 'C', ml, n, k, Aii, lda, Tii, ldt, Aij, lda, work, ldwork );

	free( tau );
	return 0;
*/
