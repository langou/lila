#include "lila.h"

int lila_dgeqr2_w02a( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork ){

	int info ; 
	double *tau, *work1;
	double *Aii, *Tii, *Qii;
	int ml, lwork1;

	tau = work;
	work1 = work+n;
	lwork1 = lwork-n;

	Aii = A + i*lda + i;
	Qii = Q + i*ldq + i;
	Tii = T + i*ldt + i;

	ml = m - i;

  	info = LAPACKE_dgeqrf_work( LAPACK_COL_MAJOR, ml, n, Aii, lda, tau, work1, lwork1 );

	info = LAPACKE_dlarft_work( LAPACK_COL_MAJOR, 'F', 'C', ml, n, Aii, lda, tau, Tii, ldt);
	
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml, n, Aii, lda, Qii, ldq );

	info = LAPACKE_dorgqr_work( LAPACK_COL_MAJOR, ml, n, n, Qii, ldq, tau, work1, lwork1 );

	return 0;

}
