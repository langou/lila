#include "qr2.h"

int lapack_mod_dlarfb_lnfc_bz( int m, int n, int k, double *V, int ldv, double *T, int ldt, double *C, int ldc ){

	double *C1, *C2;
	double *V1, *V2;

	C1 = C;
	C2 = C + k;

	V1 = V;
	V2 = V + k;

//	C1 = C1 + V2^T * C2
	cblas_dgemm( CblasColMajor, CblasTrans, CblasNoTrans, k, n, m-k, (+1.0e+00), V2, ldv, C2, ldc, (+0.0e+00), C1, ldc );

//	C1 = T * C1
	cblas_dtrmm( CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, k, n, (+1.0e+00), T, ldt, C1, ldc );

//	C2 = C2 - V2 * C1
	cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, m-k, n, k, (-1.0e+00), V2, ldv, C1, ldc, (+1.0e+00), C2, ldc );

//	C1 = - V1 * C1
	cblas_dtrmm( CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, k, n, (-1.0e+00), V1, ldv, C1, ldc );

	return 0;

}
