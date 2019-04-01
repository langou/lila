#include "qr3.h"

int qr3_test_qq_orth_1( double *norm_orth_1, int m, int n, double *Q, int ldq  ){

	double *work;
	int info, lwork;

	lwork = n;
	work  = (double *) malloc( lwork * n * sizeof(double));
	info  = LAPACKE_dlaset( LAPACK_COL_MAJOR, 'A', n, n, (0e+00), (1e+00), work, lwork );
	cblas_dsyrk( CblasColMajor, CblasUpper, CblasTrans, n, m, 1.0e+00, Q, ldq, -1.0e+00, work, lwork );
	(*norm_orth_1) = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', n, n, work, lwork, NULL );
	free( work );

	return 0;
}
