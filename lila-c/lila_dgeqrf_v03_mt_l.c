#include "lila.h"

int lila_dgeqrf_v03_mt_l( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *work, int lwork ){

	double *tau=NULL, *Aii;
	int ml, info; 

	tau = work + n;

	Aii = A + i*lda + i;

	ml = m - i;

  	info = LAPACKE_dgeqrf_work( LAPACK_COL_MAJOR, ml, n, Aii, lda, tau, work, lwork );
	info = lila_dlarft_w03( m, n, i, mt, A, lda, T, ldt, tau);

	return 0;

}
