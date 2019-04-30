#include "qr3.h"

int qr3_larft( int m, int n, double *A, int lda, double *T, int ldt, double *tau ){

	int n1, n2, i, j;
	double *A11, *A21, *A22, *A32, *A31;
	double *T11, *T12, *T22;
	double *tau1, *tau2;

	if ( n <= 1){

		*T = *tau;

	} else {

	n1 = n/2;
	n2 = n-n1;

	A11 = A;
	A21 = A+n1;
	A31 = A+n;
	A22 = A+n1*lda+n1;
	A32 = A+n1*lda+n;

	T11 = T;
	T12 = T+n1*ldt;
	T22 = T+n1*ldt+n1;

	tau1 = tau;
	tau2 = tau+n1;

	qr3_larft( m, n1, A11, lda, T11, ldt, tau1 );
	qr3_larft( m-n1, n2, A22, lda, T22, ldt, tau2 );

	for (i=0;i<n1;i++) for (j=0;j<n2;j++) T12[i+j*ldt] = A21[j+i*lda];
	cblas_dtrmm ( CblasColMajor, CblasRight, CblasLower, CblasNoTrans, CblasUnit, n1, n2, 1.0e+00, A22, lda, T12, ldt); 
	cblas_dgemm( CblasColMajor, CblasTrans, CblasNoTrans, n1, n2, m-n, +1.0e+00, A31, lda, A32,lda, +1.0e+00, T12, ldt);
	cblas_dtrmm ( CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, n1, n2, -1.0e+00, T, ldt, T12, ldt); 
	cblas_dtrmm ( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, n1, n2, +1.0e+00, T22, ldt, T12, ldt); 

	}	

	return 0;

}
