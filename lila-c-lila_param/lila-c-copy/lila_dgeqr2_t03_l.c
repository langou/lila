#include "lila.h"

int lila_dgeqr2_t03_l( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *tau, double *work, int lwork ){

	double *Aii, *Tki;
	int k, ml, info; 

	ml = m-i;
	k  = i%mt;

	Aii = A + i + i*lda;
	Tki = T + k + i*ldt;
	tau = tau + i;

	info = LAPACKE_dlarft_work( LAPACK_COL_MAJOR, 'F', 'C', ml, n, Aii, lda, tau, Tki, ldt);

	return 0;

}
