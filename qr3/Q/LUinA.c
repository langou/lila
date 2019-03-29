#include "lila.h"

int LUinA( int n, double *L, int ldl, double *A, int lda ){

	int info; 
	double *work;

	work  = (double *) malloc( n * n * sizeof(double));

	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', n, n, A, lda, work, n );
	info = LAPACKE_dlaset( LAPACK_COL_MAJOR, 'U', n, n, (0.0e+00), (1.0e+00), work, n );

	cblas_dtrmm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, n, n, (-1.0e+00), L, ldl, work, n );

	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', n, n, work, n, A, lda );
	free( work );

	return 0;

}


