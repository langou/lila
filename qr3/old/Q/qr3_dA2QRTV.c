#include "qr3.h"

int qr3_dA2QRTV( int m, int n, double *A, int lda, double *Q, int ldq, double *T, int ldt ){

	double *A11, *A12, *A21, *A22;
	double *Q11, *Q12, *Q21, *Q22;
	double *T11, *T12, *T22;

	int n1, n2, i, j;

	if( n <= 1 ){

		LAPACKE_dlarfg_work( m, A, A+1, 1, T);
		for( i=1; i<m; i++)  Q[i] = A[i] * (-(*T));
		(*Q) = (+1.0e+00) - (*T);

	} else {

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

		T11 = T;
		T12 = T+n1*ldt;
		T22 = T+n1+n1*ldt;

		qr3_dA2QRTV( m, n1, A11, lda, Q11, ldq, T11, ldt );

		// push
		for (i=0;i<n1;i++) for (j=0;j<n2;j++) T12[i+j*ldt] = A12[i+j*lda];
		cblas_dtrmm( CblasColMajor, CblasLeft, CblasLower, CblasTrans, CblasUnit, n1, n2, (+1.0e+00), A11, lda, T12, ldt );
		cblas_dgemm( CblasColMajor, CblasTrans, CblasNoTrans, n1, n2, m-n1, (+1.0e+00), A21, lda, A22, lda, (+1.0e+00), T12, ldt );
		cblas_dtrmm( CblasColMajor, CblasLeft, CblasUpper, CblasTrans, CblasNonUnit, n1, n2, (+1.0e+00), T11, ldt, T12, ldt );
		cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, m-n1, n2, n1, (-1.0e+00), A21, lda, T12, ldt, (+1.0e+00), A22, lda );
		cblas_dtrmm( CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, n1, n2, (+1.0e+00), A11, lda, T12, ldt );
		for (i=0;i<n1;i++) for (j=0;j<n2;j++) A12[i+j*lda] -= T12[i+j*ldt];

		qr3_dA2QRTV( m-n1, n2, A22, lda, Q22, ldq, T22, ldt );

		// pop
		cblas_dgemm( CblasColMajor, CblasTrans, CblasNoTrans, n1, n2, m-n1, (+1.0e+00), A21, lda, Q22, ldq, (+0.0e+00), Q12, ldq );
		cblas_dtrmm( CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, n1, n2, (+1.0e+00), T11, ldt, Q12, ldq );
		cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, m-n1, n2, n1, (-1.0e+00), A21, lda, Q12, ldq, (+1.0e+00), Q22, ldq );
		cblas_dtrmm( CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, n1, n2, (-1.0e+00), A11, lda, Q12, ldq );

		// connect
		for (i=0;i<n1;i++) for (j=0;j<n2;j++) T12[i+j*ldt] = A21[j+i*lda];
		cblas_dtrmm ( CblasColMajor, CblasRight, CblasLower, CblasNoTrans, CblasUnit, n1, n2, 1.0e+00, A22, lda, T12, ldt); 
		cblas_dgemm ( CblasColMajor, CblasTrans, CblasNoTrans, n1, n2, m-n, +1.0e+00, A21+n2, lda, A22+n2,lda, +1.0e+00, T12, ldt);
		cblas_dtrmm ( CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, n1, n2, -1.0e+00, T11, ldt, T12, ldt); 
		cblas_dtrmm ( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, n1, n2, +1.0e+00, T22, ldt, T12, ldt); 

	}

	return 0;

}
