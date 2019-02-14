#include "lila.h"

int lila_wsq_dgeqrf_w03_mt_l( int panel, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork ){

	double *tau=NULL, *Aii, *Qii;
	int ml, info; 

	tau = work + n;

	Aii = A + i*lda + i;
	Qii = Q + i*ldq + i;

	ml = m - i;

  	info = LAPACKE_dgeqrf_work( LAPACK_COL_MAJOR, ml, n, Aii, lda, tau, work, lwork );
//	info = lila_dlarft_w03( m, n, i, mt, A, lda, T, ldt, tau);
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml, n, Aii, lda, Qii, ldq );
	info = LAPACKE_dorgqr_work( LAPACK_COL_MAJOR, ml, n, n, Qii, ldq, tau, work, lwork );


	return n;

}
