#include "lila.h"

int lila_qr_test( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *As, double normA ){
// int lila_test_qr_repres_1( int m, int n, int i, double *A, int lda, double *Q, int ldq, double *R, int ldr, double normA ){

	double *Qii, *Aii, *As_ii;
	double norm_repres, *work;
	int ml, info, lwork, jj, ii;

	ml = m-i;

	As_ii = As+i+i*lda;
	Aii   =  A+i+i*lda;

	lwork = ml*n;
	work  = (double *) malloc(ml * n * sizeof(double));
	info  = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml, n, Qii, ldq, work, ml );
	cblas_dtrmm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, ml, n, (1.0e+00), Aii, lda, work, ml );
 	for(ii = 0; ii < ml; ii++) for(jj = 0; jj < n; jj++) work[ ii+jj*ml ] -= As_ii[ ii+jj*lda ];
	norm_repres = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', ml, n, work, ml, NULL );
	norm_repres = norm_repres / normA;
	free( work );

	return norm_repres;
}

