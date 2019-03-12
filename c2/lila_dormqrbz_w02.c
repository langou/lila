#include "lila.h"

int lila_dormqrbz_w02( int m, int n, int k, int i, int j, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork ){

	double *Aii, *Qij;
	int ml;
	double *Tii;
	int ldwork;
	
	Aii = A + i + i*lda;
	Qij = Q + i + j*ldq;
	Tii = T + i + i*ldt;

	ml = m - i;
	ldwork = k;


	printf("ldq = %d, k = %d\n",ldq,k);
 	cblas_dgemm( CblasColMajor, CblasTrans, CblasNoTrans, k, n, m-k-i, (1.0e+00), Aii+k, lda, Qij+k, lda, (0.0e+00), work, ldwork );
	cblas_dtrmm( CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, k, n, (1.0e+00), Tii, ldt, work, ldwork );
	LAPACKE_dlacpy( LAPACK_COL_MAJOR, 'A', k, n, work, ldwork, Qij, ldq );
	cblas_dtrmm( CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, k, n, (-1.0e+00), Aii, lda, Qij, ldq );
	cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, m-k-i, n, k, (-1.0e+00), Aii+k, lda, work, ldwork, (1.0e+00), Qij+k, ldq );// not sure why we are running into an issue on param 13

	return 0;

}
