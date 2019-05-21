#include "A2VR.h"

int dgeqr3R_UT( int m, int n, double *A, int lda, double *T, int ldt, double *R, int ldr ){

	int n1, n2, i, j, info;
	double *A11, *A12, *A21, *A22;
	double *T11, *T12, *T22;
	double *R11, *R12, *R22;

	if ( n == 1){

		double tau;
		info = LAPACKE_dlarfg_work( m, (&A[0]), &(A[1]), 1, &tau);
		(*R) = (*A); 
		(*T) = 1/tau;

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

		R11 = R;
		R12 = R+n1*ldr;
		R22 = R+n1*(1+ldr);

	dgeqr3R_UT( m, n1, A11, lda, T11, ldt, R11, ldr );

		// Using T12 as a workspace
		if( A11 == R11 ){

		for (i=0;i<n1;i++) for (j=0;j<n2;j++) T12[i+j*ldt] = A12[j*lda+i];
		cblas_dtrmm ( CblasColMajor, CblasLeft, CblasLower, CblasTrans, CblasUnit, n1, n2, 1.0e+00, A11, lda, T12, ldt); 
		cblas_dgemm( CblasColMajor, CblasTrans, CblasNoTrans, n1, n2, m-n1, +1.0e+00, A21, lda, A22, lda, +1.0e+00, T12, ldt);
		cblas_dtrsm ( CblasColMajor, CblasLeft, CblasUpper, CblasTrans, CblasNonUnit, n1, n2, 1.0e+00, T11, ldt, T12, ldt); 
		cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, m-n1, n2, n1, -1.0e+00, A21, lda, T12, ldt, +1.0e+00, A22, lda);
		cblas_dtrmm ( CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, n1, n2, 1.0e+00, A11, lda, T12, ldt); 
		for (j=0;j<n2;j++) for (i=0;i<n1;i++) R12[j*ldr+i] = A12[j*lda+i] - T12[j*ldt+i]; 

		// Using R12 as a workspace
		} else {              // ( A11 == T11 ) || ( A11 != T11 & A11 != R11 ) // 

		for (i=0;i<n1;i++) for (j=0;j<n2;j++) R12[i+j*ldr] = A12[j*lda+i];
		cblas_dtrmm ( CblasColMajor, CblasLeft, CblasLower, CblasTrans, CblasUnit, n1, n2, 1.0e+00, A11, lda, R12, ldr); 
		cblas_dgemm( CblasColMajor, CblasTrans, CblasNoTrans, n1, n2, m-n1, +1.0e+00, A21, lda, A22, lda, +1.0e+00, R12, ldr);
		cblas_dtrsm ( CblasColMajor, CblasLeft, CblasUpper, CblasTrans, CblasNonUnit, n1, n2, 1.0e+00, T11, ldt, R12, ldr); 
		cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, m-n1, n2, n1, -1.0e+00, A21, lda, R12, ldr, +1.0e+00, A22, lda);
		cblas_dtrmm ( CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, n1, n2, 1.0e+00, A11, lda, R12, ldr); 
		for (j=0;j<n2;j++) for (i=0;i<n1;i++) R12[j*ldr+i] = A12[j*lda+i] - R12[j*ldr+i]; 

		}

	dgeqr3R_UT( m-n1, n2, A22, lda, T22, ldt, R22, ldr );

		for (i=0;i<n1;i++) for (j=0;j<n2;j++) T12[i+j*ldt] = A21[j+i*lda];
		cblas_dtrmm ( CblasColMajor, CblasRight, CblasLower, CblasNoTrans, CblasUnit, n1, n2, 1.0e+00, A22, lda, T12, ldt); 
		cblas_dgemm( CblasColMajor, CblasTrans, CblasNoTrans, n1, n2, m-n, +1.0e+00, A21+n2, lda, A22+n2,lda, +1.0e+00, T12, ldt);

	}


	return 0;

}
