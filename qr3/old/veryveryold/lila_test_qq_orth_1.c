#include "lila.h"

double lila_test_qq_orth_1( int m, int n, int i, double *Q, int ldq  ){

	double *Qii, norm_orth_1, *work;
	int ml, info, lwork;

	ml = m - i;
	Qii = Q + i + i*ldq;

	lwork = n;
	work  = (double *) malloc( lwork * n * sizeof(double));
	info  = LAPACKE_dlaset( LAPACK_COL_MAJOR, 'A', n, n, (0e+00), (1e+00), work, lwork );
	cblas_dsyrk( CblasColMajor, CblasUpper, CblasTrans, n, ml, 1.0e+00, Qii, ldq, -1.0e+00, work, lwork );
	norm_orth_1 = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', n, n, work, lwork, NULL );
	free( work );

	return norm_orth_1;
}

