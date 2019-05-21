#include "V2T.h"

int V2T_test_qr_repres_1( double *norm_repres, int m, int n, double *A, int lda, double *Q, int ldq, double *R, int ldr ){

	double normA, *work;
	int ii, jj;

	normA = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', m, n, A, lda, NULL );

	work  = (double *) malloc(m * n * sizeof(double));
	LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, n, Q, ldq, work, m );
	cblas_dtrmm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, m, n, (1.0e+00), R, ldr, work, m );
 	for(ii = 0; ii < m; ii++) for(jj = 0; jj < n; jj++) work[ ii+jj*m ] -= A[ ii+jj*lda ];
	(*norm_repres) = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', m, n, work, m, NULL );
	free( work );

	(*norm_repres) = (*norm_repres) / normA;

	return 0;
}

