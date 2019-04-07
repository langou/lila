#include "qr3.h"

int dN2T( int n, double *tau, double *T, int ldt ){

	int n1, n2; 
	double *T11, *T12, *T22;
	double *tau1, *tau2;

	if ( n <= 1 ) {

// those three lines below are correct, we keep all three for now
// the two last formulae (based on 2/T ) are significantly less accurate than the formula based on tau

		(*T) = (*tau);
//		if ( (*T)!= 0.0e+00 ) (*T) = (+2.0e+00) / (*T); 
//		(*T) = ( (*T)== 0.0e+00 ? 0.0e+00 : (+2.0e+00) / (*T) );

	} else {

		n1 = n / 2;
		n2 = n - n1;

		T11 = T;
		T12 = T+n1*ldt;
		T22 = T+n1+n1*ldt;

		tau1 = tau;
		tau2 = tau+n1;

		dN2T( n1, tau1, T11, ldt );
		dN2T( n2, tau2, T22, ldt );

		cblas_dtrmm( CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit,  n1, n2, (-1.0e+00), T11, ldt, T12, ldt );
		cblas_dtrmm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, n1, n2, (+1.0e+00), T22, ldt, T12, ldt );

	}

	return 0;

}
