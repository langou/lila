#include "lila.h"

int lila_wsq_dgeqrf_w03_mt_l( int panel, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork ){

//	tau = work + n;
//  	info = LAPACKE_dgeqrf_work( LAPACK_COL_MAJOR, ml, n, Aii, lda, tau, work, lwork );
//	info = lila_dlarft_w03( m, n, i, mt, A, lda, T, ldt, tau);
//	info = LAPACKE_dorgqr_work( LAPACK_COL_MAJOR, ml, n, n, Qii, ldq, tau, work, lwork );

	return n; // or is it n*n

}
