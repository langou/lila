#include "check.h"

int check_qq_orth( double *norm_orth_1, int m, int n, double *Q, int ldq  ){

	double *work;

	work  = (double *) malloc( n * n * sizeof(double));
	LAPACKE_dlaset( LAPACK_COL_MAJOR, 'A', n, n, (0e+00), (1e+00), work, n );
	cblas_dsyrk( CblasColMajor, CblasUpper, CblasTrans, n, m, 1.0e+00, Q, ldq, -1.0e+00, work, n );
	(*norm_orth_1) = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', n, n, work, n, NULL );
	free( work );

	return 0;
}
