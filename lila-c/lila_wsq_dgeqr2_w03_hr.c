#include "lila.h"

int lila_wsq_dgeqr2_w03_hr( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork ){

	double *Aii, *Qii, *Tki;
	double *S;
	int info, k, ml;

//	int *S;
//	S = (int *) malloc(n * sizeof(int));
	S = work;

	ml = m - i;
	k  =  i%mt;

	Aii = A + i + i*lda;
	Qii = Q + i + i*ldq;
	Tki = T + k + i*ldt;
	
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml, n, Aii, lda, Qii, ldq ); 
	cblas_dsyrk( CblasColMajor, CblasUpper, CblasTrans, n, ml, 1.0e+00, Qii, ldq, 0e+00, Aii, lda );
	info = LAPACKE_dpotrf( LAPACK_COL_MAJOR, 'U', n, Aii, lda ); 
	cblas_dtrsm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, ml, n, 1.0e+00, Aii, lda, Qii, ldq );

	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml-n, n, Qii+n, ldq, Aii+n, lda ); 
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', n, n, Qii, ldq, Tki, ldt ); 
	lila_dorgh2_3( ml, n, Aii, lda, Tki, ldt, Qii, ldq, work, lwork, S );

//	free(S);

	return 0;

}
