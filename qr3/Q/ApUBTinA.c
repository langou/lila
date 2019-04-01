#include "qr3.h"

//	A is m-by-n
//	U is m-by-m and upper triangular
//	B is n-by-m
//
//	ApUBTinA computes A+UB^T in-place of A taking into account the fact that U is upper triangular
//
//	algorithm is recursive and breaks down A, B, and U in four submatrices
//
//	A11 is m1-by-n1
//	A12 is m1-by-n2
//	A21 is m2-by-n1
//	A22 is m2-by-n2
//
//	U11 is m1-by-m1
//	U12 is m1-by-m2
//	U22 is m2-by-m2
//
//	B11 is n1-by-m1
//	B12 is n1-by-m2
//	B21 is n2-by-m1
//	B22 is n2-by-m2
//
//	we want to perform
//		A = A + U B^T
//
//	we want to perform
//		( A11 A12 ) = ( A11 A12 ) + ( U11 U12 ) ( B11^T B21^T )
//		( A21 A22 )   ( A21 A22 )   ( 0   U22 ) ( B12^T B22^T )
//
//	we want to perform
//		A11 = A11 + U11 * B11^T + U12 * B12^T
//		A12 = A12 + U11 * B21^T + U12 * B22^T
//		A21 = A21 + U22 * B12^T (myself)
//		A22 = A22 + U22 * B22^T (myself)
//
//	we want to perform
//		A11 = A11 + U11 * B11^T (myself)
//		A11 = A11 + U12 * B12^T (gemm)
//		A12 = A12 + U11 * B21^T (myself)
//		A12 = A12 + U12 * B22^T (gemm)
//		A21 = A21 + U22 * B12^T (myself)
//		A22 = A22 + U22 * B22^T (myself)

int ApUBTinA( int m, int n, double *A, int lda, double *U, int ldu, double *B, int ldb ){

	int m1, m2, n1, n2;
	double *A11, *A12, *A21, *A22;
	double *B11, *B12, *B21, *B22;
	double *U11, *U12, *U22;

	if (( m <= 1 )&&( n<=1 )){

		(*A) = (*A) + (*U) * (*B);

	} else {


		m1 = m/2;
		m2 = m-m1;

		n1 = n/2; 
		n2 = n-n1;

		A11 = A;
		A12 = A + n1*lda;
		A21 = A + m1;
		A22 = A + m1 + n1*lda;

		B11 = B;
		B12 = B + m1*ldb;
		B21 = B + n1;
		B22 = B + n1 + m1*ldb;

		U11 = U;
		U12 = U + m1*ldu;
		U22 = U + m1 + m1*ldu;

		if(( m1 > 0 )&&( n1 > 0 )) ApUBTinA( m1, n1, A11, lda, U11, ldu, B11, ldb );

		if(( m1 > 0 )&&( n1 > 0 )&&( m2 > 0 )) cblas_dgemm( CblasColMajor, CblasNoTrans, CblasTrans, m1, n1, m2, (+1.0e+00), U12, ldu, B12, ldb, (+1.0e+00), A11, lda );

		if(( m1 > 0 )&&( n2 > 0 )) ApUBTinA( m1, n2, A12, lda, U11, ldu, B21, ldb );

		if(( m1 > 0 )&&( n2 > 0 )&&( m2 > 0 )) cblas_dgemm( CblasColMajor, CblasNoTrans, CblasTrans, m1, n2, m2, (+1.0e+00), U12, ldu, B22, ldb, (+1.0e+00), A12, lda );

		if(( m2 > 0 )&&( n1 > 0 )) ApUBTinA( m2, n1, A21, lda, U22, ldu, B12, ldb );

		if(( m2 > 0 )&&( n2 > 0 )) ApUBTinA( m2, n2, A22, lda, U22, ldu, B22, ldb );

	}

	return 0;

}
