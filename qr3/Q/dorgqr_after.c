#include "qr3.h"

int dorgqr_after( int m, int n, int k, double *A, int lda, double *T, int ldt, double *Q, int ldq ){

	double *A1, *A2, *A3;
	double *Q1, *Q2, *Q3;
	int i, j;

	A1 = A;
	A2 = A + k;
	A3 = A + n;

	Q1 = Q;
	Q2 = Q + k;
	Q3 = Q + n;

	// copy V2^T into Q1
//	LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', k, n-k, A2, lda, Q1, ldq );
	for ( i = 0; i < k; i++ ) for ( j = 0; j < n-k; j++ ) Q1[i+j*ldq] = A2[j+i*lda];

	// - T * V2^T -- note putting the minus sign on this term
	cblas_dtrmm( CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit,  k, n-k, (-1.0e+00), T, ldt, Q1, ldq );

	// V2 * ( T * V2^T )
	cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, n-k, n-k, k, (+1.0e+00), A2, lda, Q1, ldq, (+0.0e+00), Q2, ldq);

	// Adding the identity block to Q2
	for( j = 0; j < n-k; j++ ){ Q2[ j + j*ldq ] += (+1.0e+00); }

	// V3 * ( T * V2^T )
	cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, m-n, n-k, k, (+1.0e+00), A3, lda, Q1, ldq, (+0.0e+00), Q3, ldq);

	// V1 * ( T * V2^T ) 
	cblas_dtrmm( CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, k, n-k, (+1.0e+00), A1, lda, Q1, ldq );

	return 0;

}
