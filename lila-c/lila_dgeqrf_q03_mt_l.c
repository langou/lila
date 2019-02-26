#include "lila.h"

int lila_dgeqrf_q03_mt_l( int panel, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork ){

	double *tau=NULL, *Aii, *Qii;
	int ml, info, lwork1; 

	ml     = m - i;
	tau    = work;
	work   = work + n;
	lwork1 = lwork-n;

	Aii = A + i*lda + i;
	Qii = Q + i*ldq + i;

	info = lila_dT2tau_w03( m, n, i, mt, T, ldt, tau);
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml, n, Aii, lda, Qii, ldq );
	info = LAPACKE_dorgqr_work( LAPACK_COL_MAJOR, ml, n, n, Qii, ldq, tau, work, lwork1 );


	return 0;

}
