#include "check.h"

int check_q2A_repres( double *norm_repres, int m, int n, int k, double *A, int lda, double *Q2, int ldq2 ){

	double normA, *work;

	normA = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', m, k, A, lda, NULL );

	work  = (double *) malloc( (n-k) * k * sizeof(double));
	cblas_dgemm( CblasColMajor, CblasTrans, CblasNoTrans, n-k, k, m, (+1.0e+00), Q2, ldq2, A, lda, (+0.0e+00), work, n-k );
	(*norm_repres) = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', n-k, k, work, n-k, NULL );

	(*norm_repres) = (*norm_repres) / normA;

	free( work );

	return 0;
}

