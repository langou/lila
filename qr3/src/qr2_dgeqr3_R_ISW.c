#include "qr2.h"

int qr2_dgeqr3_R_ISW( int m, int n, double *A, int lda, double *T, int ldt, double *R, int ldr ){

	int n1, n2, i, j, ldw, info;
	double *A11, *A12, *A21, *A22;
	double *T11, *T12, *T22;
	double *R11, *R12, *R22;
	double *W;

	if ( n == 1){

//		double tau;
//		info = LAPACKE_dlarfg_work( m, A, A+1, 1, &tau);
//		(*R) = (*A); 
//		(*T) = tau;

		(*R) = (*A); 
		info = LAPACKE_dlarfg_work( m, R, A+1, 1, T);

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

	qr2_dgeqr3_R( m, n1, A11, lda, T11, ldt, R11, ldr );

		if( A11 == R11 ){ 
			W = T12; ldw = ldt; 
		} else { 
			W = R12; ldw = ldr; 
		}

		for (i=0;i<n1;i++) for (j=0;j<n2;j++) W[i+j*ldw] = A12[j*lda+i];
		cblas_dtrmm ( CblasColMajor, CblasLeft, CblasLower, CblasTrans, CblasUnit, n1, n2, 1.0e+00, A11, lda, W, ldw); 
		cblas_dgemm( CblasColMajor, CblasTrans, CblasNoTrans, n1, n2, m-n1, +1.0e+00, A21, lda, A22, lda, +1.0e+00, W, ldw);
		cblas_dtrmm ( CblasColMajor, CblasLeft, CblasUpper, CblasTrans, CblasNonUnit, n1, n2, 1.0e+00, T11, ldt, W, ldw); 
		cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, m-n1, n2, n1, -1.0e+00, A21, lda, W, ldw, +1.0e+00, A22, lda);
		cblas_dtrmm ( CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, n1, n2, 1.0e+00, A11, lda, W, ldw); 
		for (j=0;j<n2;j++) for (i=0;i<n1;i++) R12[j*ldr+i] = A12[j*lda+i] - W[j*ldw+i]; 

	qr2_dgeqr3_R_ISW( m-n1, n2, A22, lda, T22, ldt, R22, ldr );
 

	}


	return 0;

}

