#include "qr2.h"

int lapack_ref_dlarfb_lnfc( int m, int n, int k, double *V, int ldv, double *T, int ldt, double *C, int ldc, double *W ){

	int i, j;
	double *C1, *C2;
	double *V1, *V2;

	C1 = C;
	C2 = C + k;

	V1 = V;
	V2 = V + k;

//	W = C1
	LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', k, n, C1, ldc, W, k );

//	W = V1^T * W
	cblas_dtrmm( CblasColMajor, CblasLeft, CblasLower, CblasTrans, CblasUnit, k, n, (+1.0e+00), V1, ldv, W, k );

//	W = W + V2^T * C2
	cblas_dgemm( CblasColMajor, CblasTrans, CblasNoTrans, k, n, m-k, (+1.0e+00), V2, ldv, C2, ldc, (+1.0e+00), W, k );

//	W = T * W
	cblas_dtrmm( CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, k, n, (+1.0e+00), T, ldt, W, k );

//	C2 = C2 - V2 * W
	cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, m-k, n, k, (-1.0e+00), V2, ldv, W, k, (+1.0e+00), C2, ldc );

//	W = - V1 * W
	cblas_dtrmm( CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, k, n, (-1.0e+00), V1, ldv, W, k );

//	C1 = C1 + W
	for(j = 0; j < n; j++) for(i = 0; i < k; i++) C1[ i+j*ldc ] += W[ i+j*k ];

	return 0;

}
