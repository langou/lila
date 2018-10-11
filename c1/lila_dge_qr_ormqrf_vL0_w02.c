#include "lila.h"

int lila_dge_qr_ormqrf_vL0( int m, int n, int k, int i, int j, int mt, double *A, int lda, double *T, int ldt, double *work, int lwork ){

	double *tau;
	int info, jj;
	double *Aii, *Tii, *Aij;
	int ml;
	double *V;
	double normV_square;
	int vl;
	
	Aii = A + i*lda + i;
	Aij = A + j*lda + i;
	Tii = T + i*ldt + i;
	ml = m - i;

	tau = (double *) malloc( k * sizeof(double));
//	for(jj = 0; jj < k; jj++) tau[jj] = T[jj+jj*ldt];

	for(jj = 0, vl=ml-1, V = Aii+1; jj < k; jj++,vl--,V+=(lda+1) ){
		normV_square = ( 1.0e+00 ) + cblas_ddot( vl, V, 1, V, 1 );
		tau[jj] = ( 2.0e+00 ) / normV_square ;
	}

	info = LAPACKE_dormqr_work( LAPACK_COL_MAJOR, 'L', 'T', ml, n, k, Aii, lda, tau, Aij, lda, work, lwork );

	free( tau );

	return 0;

}
