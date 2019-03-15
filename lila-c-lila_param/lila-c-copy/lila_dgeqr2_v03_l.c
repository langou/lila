#include "lila.h"

int lila_dgeqr2_v03_l( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *work, int lwork ){

	double *tau=NULL, *Aii, *Tki;
	int k, ml, info, lwork1; 

	tau = work;
	work = work + n;
	lwork1 = lwork-n;

	ml = m-i;
	k  = i%mt;

	Aii = A + i + i*lda;
	Tki = T + k + i*ldt;

  	info = LAPACKE_dgeqrf_work( LAPACK_COL_MAJOR, ml, n, Aii, lda, tau, work, lwork1 ); 
	info = LAPACKE_dlarft_work( LAPACK_COL_MAJOR, 'F', 'C', ml, n, Aii, lda, tau, Tki, ldt);

	return 0;

}
