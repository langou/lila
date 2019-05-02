#include "qr3.h"

int qr3_dA2QRTV( int m, int n, double *A, int lda, double *Q, int ldq, double *R, int ldr, double *T, int ldt, double *V, int ldv ){

	double A1, A2;
	double Q11, Q22, Q12, Q21;
	double R1, R2;
	double T1, T2;
	double V1, V2;

	int n1, n2;

	n1 = n / 2;
	n2 = n - n1;

	A11 = A;
	A21 = A+n1;
	A12 = A+n1*lda;
	A22 = A+n1+n1*lda;

	Q11 = Q;
	Q21 = Q+n1;
	Q12 = Q+n1*ldq;
	Q22 = Q+n1+n1*ldq;

	R11 = R;
	R21 = R+n1;
	R12 = R+n1*ldr;
	R22 = R+n1+n1*ldr;

	T11 = T;
	T21 = T+n1;
	T12 = T+n1*ldt;
	T22 = T+n1+n1*ldt;

	V11 = V;
	V21 = V+n1;
	V12 = V+n1*ldv;
	V22 = V+n1+n1*ldv;

	qr3_dA2QRTV( m, n1, A, lda, Q, ldq, R, ldr, T, ldt, V, ldv );

	cblas_dgemm( CblasColMajor, CblasTrans, CblasNoTrans, n1, n2, m-n1, (+1.0e+00), Q21, ldq, Q22, ldq, (+0.0e+00), Q12, ldq );
	cblas_dtrmm( CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, n1, n2, (+1.0e+00), T11, ldt, Q12, ldq );
	cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, m-n1, n2, n1, (-1.0e+00), Q21, ldq, Q12, ldq, (+1.0e+00), Q22, ldq );
	cblas_dtrmm( CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, n1, n2, (-1.0e+00), Q11, ldq, Q12, ldq );

	qr3_dA2QRTV( m, n1, A, lda, Q, ldq, R, ldr, T, ldt, V, ldv );




	return 0;

}
