#include "lila.h"

int lila_dgeqr2_w03b( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork ){

	int info ; 
	double *Aii, *Ajj, *Qii;
	int ml, i1, j1;

	Aii = A + i*lda + i;
	Qii = Q + i*ldq + i;
	
	ml = m - i;

//  This block computes the Cholesky QR
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml, n, Aii, lda, Qii, ldq ); 

	cblas_dsyrk( CblasColMajor, CblasUpper, CblasTrans, n, ml, 1.0e+00, Qii, ldq, 0e+00, Aii, lda );
	info = LAPACKE_dpotrf( LAPACK_COL_MAJOR, 'U', n, Aii, lda ); 
	cblas_dtrsm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, ml, n, 1.0e+00, Aii, lda, Qii, ldq );

//   LU and construct T

//	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml, n, Qii, ldq, Aii, lda );
//	lila_dgetrf_b( m, n, i, mt, A, lda, Q, ldq, R, lwork );
//	info = lila_dlarft_w03_b( m, n, i, mt, A, lda, T, ldt );

	lila_dorghr( m, n, i, mt, A, lda, T, ldt, Q, ldq, work, lwork );

	return 0;

}
