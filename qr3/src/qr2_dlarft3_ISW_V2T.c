#include "qr2.h"

int qr2_dlarft3_ISW_V2T( int m, int n, double *A, int lda, double *T, int ldt, double *tau ){

	int n1, n2;
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

	{ int i,j; for(i=0;i<n1;i++){ for(j=0;j<i;j++){ T11[j+i*ldt] = A11[i+j*lda];}} }
	qr2_aux_dV2N( n1, T11, ldt );
	cblas_dsyrk( CblasColMajor, CblasUpper, CblasTrans, n1, m-n1, (+1.0e+00), A11+n1, lda, (+1.0e+00), T11, ldt );
	qr2_aux_dN2T( n1, tau1, T11, ldt );

	qr2_dlarft3_ISW_V2T( m-n1, n2, A22, lda, T22, ldt, tau2 );

	}	

	return 0;

}
