#include "A2Q.h"

int dgeqr3_Q( int m, int n, double *A, int lda, double *T, int ldt ){

	int n1, n2, i, j, info;
	double *A11, *A12, *A21, *A22;
	double *T11, *T12, *T22;

	if ( n == 1){

		double tau;
		info = LAPACKE_dlarfg_work( m, A, A+1, 1, &tau);
		(*T) = tau;

	} else {

		n1 = n/2;
		n2 = n-n1;

		A11 = A;
		A12 = A+n1*lda;
		A21 = A+n1;
		A22 = A+n1*(1+lda);

		T11 = T;
		T12 = T+n1*ldt;
		T22 = T+n1*(1+ldt);

	dgeqr3_Q( m, n1, A11, lda, T11, ldt );

		for (i=0;i<n1;i++) for (j=0;j<n2;j++) T12[i+j*ldt] = A12[j*lda+i];
		cblas_dtrmm ( CblasColMajor, CblasLeft, CblasLower, CblasTrans, CblasUnit, n1, n2, 1.0e+00, A11, lda, T12, ldt); 
		cblas_dgemm( CblasColMajor, CblasTrans, CblasNoTrans, n1, n2, m-n1, +1.0e+00, A21, lda, A22, lda, +1.0e+00, T12, ldt);
		cblas_dtrmm ( CblasColMajor, CblasLeft, CblasUpper, CblasTrans, CblasNonUnit, n1, n2, 1.0e+00, T11, ldt, T12, ldt); 
		cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, m-n1, n2, n1, -1.0e+00, A21, lda, T12, ldt, +1.0e+00, A22, lda);
	//	cblas_dtrmm ( CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, n1, n2, 1.0e+00, A11, lda, T12, ldt); 
	//	for (j=0;j<n2;j++) for (i=0;i<n1;i++) R12[j*ldr+i] = A12[j*lda+i] - T12[j*ldt+i]; 


	dgeqr3_Q( m-n1, n2, A22, lda, T22, ldt );

		for (i=0;i<n1;i++) for (j=0;j<n2;j++) T12[i+j*ldt] = A21[j+i*lda];
		cblas_dtrmm ( CblasColMajor, CblasRight, CblasLower, CblasNoTrans, CblasUnit, n1, n2, 1.0e+00, A22, lda, T12, ldt); 
		cblas_dgemm( CblasColMajor, CblasTrans, CblasNoTrans, n1, n2, m-n, +1.0e+00, A21+n2, lda, A22+n2,lda, +1.0e+00, T12, ldt);
		cblas_dtrmm ( CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, n1, n2, -1.0e+00, T11, ldt, T12, ldt); 
		cblas_dtrmm ( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, n1, n2, +1.0e+00, T22, ldt, T12, ldt); 

	}


	return 0;

}
