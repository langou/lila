#include "lila.h"

int lila_dgeqr2_w03_l( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork ){

	int info ; 
	double *tau=NULL;
	double *Aii, *Qii, *Tki;
	int k, ml;

	tau = (double *) malloc( n * sizeof(double));

	ml = m - i;
	k = i % mt;

	Aii = A + i*lda + i;
	Qii = Q + i*ldq + i;
	Tki = T + k + i*ldt;

  	info = LAPACKE_dgeqrf_work( LAPACK_COL_MAJOR, ml, n, Aii, lda, tau, work, lwork ); 
	info = LAPACKE_dlarft_work( LAPACK_COL_MAJOR, 'F', 'C', ml, n, Aii, lda, tau, Tki, ldt);
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml, n, Aii, lda, Qii, ldq ); 
	info = LAPACKE_dorgqr_work( LAPACK_COL_MAJOR, ml, n, n, Qii, ldq, tau, work, lwork );

	free( tau );

	return 0;

}
