#include "V2T.h"

int mLUinA( int n, double *A, int lda ){

	int info, n1, n2;
	double *A11, *A12, *A21, *A22;

	if ( n <= 1 ){

		(*A) = -(*A);

	} else {

		n1 = n/2; n2 = n-n1;

		A11 = A;
		A21 = A + n1; 
		A12 = A + n1*lda;
		A22 = A + n1 + n1*lda; 

//		A22 = - L22 * U22
		info = mLUinA( n2, A22, lda );

//		A22 = A22 - L21 * U12 
		cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, n2, n2, n1, (-1.0e+00), A21, lda, A12, lda, (+1.0e+00), A22, lda );

//		A12 = - L11 * U12
		cblas_dtrmm( CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, n1, n2, (-1.0e+00), A11, lda, A12, lda );

//		A21 = - L21 * U11
		cblas_dtrmm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, n2, n1, (-1.0e+00), A11, lda, A21, lda );

//		A11 = - L11 * U11
		info = mLUinA( n1, A11, lda );

	}

	return 0;

}
