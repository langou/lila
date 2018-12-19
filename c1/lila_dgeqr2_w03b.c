#include "lila.h"

int lila_dgeqr2_w03b( int m, int n, int i, int j, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork, int *S ){

	int info; 
	double *Ajj, *Qjj;
	int ml;

	Ajj = A + j*lda + j;
	Qjj = Q + j*ldq + j;
	
	ml = m - j;

//  This block computes the Cholesky QR
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml, n, Ajj, lda, Qjj, ldq ); 

	cblas_dsyrk( CblasColMajor, CblasUpper, CblasTrans, n, ml, 1.0e+00, Qjj, ldq, 0e+00, Ajj, lda );
	info = LAPACKE_dpotrf( LAPACK_COL_MAJOR, 'U', n, Ajj, lda ); 
	cblas_dtrsm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, ml, n, 1.0e+00, Ajj, lda, Qjj, ldq );

//   LU and construct T
	lila_dorghr_w03( m, n, i, j, 0, mt, A, lda, T, ldt, Q, ldq, work, lwork, S );

	return 0;

}
