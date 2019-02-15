#include "lila.h"

int lila_wsq_dgeqr2_w03_3( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork ){

//	work1 = work;	

//	info = LAPACKE_dlaset( LAPACK_COL_MAJOR, 'A', n, n, (0.0e+00), (0.0e+00), work1, n);
//	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'U', n, n, Tki, ldt, work1, n );
//	cblas_dtrmm( CblasColMajor, CblasRight, CblasLower, CblasTrans, CblasUnit, n, n, (+1.0e+00), Aii, lda, work1, n );
//	cblas_dtrmm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, ml, n, (-1.0e+00), work1, n, Qii, ldq );

	return n;

}
