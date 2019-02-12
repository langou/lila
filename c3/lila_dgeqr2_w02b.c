#include "lila.h"

int lila_dgeqr2_w02b( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork ){

	int info ; 
	double *tau=NULL;
	double *Aii, *Tii, *Qii, *R;
	int ml, lwork1;

//	tau = (double *) malloc( n * sizeof(double));
//	tau = work+n; // work+n was of the right size for w03, but w02 needs to account for i
	tau = work;
	work = work +n;
	lwork1 = lwork-n;

	Aii = A + i*lda + i;
	Qii = Q + i*ldq + i;
	Tii = T + i*ldt + i;

	R = (double *) malloc( n * n * sizeof(double));

	ml = m - i;

	cblas_dsyrk( CblasColMajor, CblasUpper, CblasTrans, n, ml, 1.0e+00, Aii, lda, 0e+00, R, n );
	info = LAPACKE_dgeqrf_work( LAPACK_COL_MAJOR, ml, n, Aii, lda, tau, work, lwork1 );
	info = LAPACKE_dlarft_work( LAPACK_COL_MAJOR, 'F', 'C', ml, n, Aii, lda, tau, Tii, ldt);
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml, n, Aii, lda, Qii, ldq );
	info = LAPACKE_dorgqr_work( LAPACK_COL_MAJOR, ml, n, n, Qii, ldq, tau, work, lwork1 );

//	free( tau );

	return 0;

}
