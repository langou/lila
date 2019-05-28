#include "qr2.h"

int qr2_larft3_ISW_V2T_UT( int m, int n, double *A, int lda, double *T, int ldt, double *tau ){

	int n1, n2;
	double *A11, *A21, *A22, *A32, *A31;
	double *T11, *T12, *T22;
	double *tau1, *tau2;

	if ( n <= 1){

		*T = +1.0e00 / *tau;

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
	qr2_dV2N( n1, T11, ldt );
	cblas_dsyrk( CblasColMajor, CblasUpper, CblasTrans, n1, m-n1, (+1.0e+00), A11+n1, lda, (+1.0e+00), T11, ldt );
	{ int i; for(i=0;i<n1;i++){ T11[i+i*ldt] = +1.0e00 / tau1[i]; } }

	qr2_larft3_ISW_V2T_UT( m-n1, n2, A22, lda, T22, ldt, tau2 );

	}	

	return 0;

}
