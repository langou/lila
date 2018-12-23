#include "lila.h"

int lila_wsq_dgeqr2_w02a( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work ){

	int info ; 
	double *tau=NULL;
	double *Aii, *Tii, *Qii;
	int ml;
	double w;
	int lwork, lwork1;

	ml = m - i;

	lwork = 0;

  	info = LAPACKE_dgeqrf_work( LAPACK_COL_MAJOR, ml, n, Aii, lda, tau, &w, -1 ); lwork1 = ((int) w); if ( lwork1 > lwork ) lwork = lwork1;

	info = LAPACKE_dorgqr_work( LAPACK_COL_MAJOR, ml, n, n, Qii, ldq, tau, &w, -1 ); lwork1 = ((int) w); if ( lwork1 > lwork ) lwork = lwork1;

	lwork += n;

	return lwork;

}
