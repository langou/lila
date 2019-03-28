#include "lila.h"

int xN2T( int n, double *tau, double *T, int ldt ){

	double *Tii, *taui;
	double *T0i;
	int info, n1, n2, ii, jj; 

	if ( n <= 1 ) {

//		(*T) = (*tau);
		(*T) = (+2.0e+00) / (*T);

	} else {

		n1 = n / 2;
		n2 = n - n1;

		Tii = T + n1 + n1*ldt;
		taui = tau + n1;

		T0i = T + n1*ldt;

		info = xN2T( n1, tau, T, ldt );
		info = xN2T( n2, taui, Tii, ldt );

		cblas_dtrmm( CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit,  n1, n2, (-1.0e+00), T, ldt, T0i, ldt );
		cblas_dtrmm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, n1, n2, (+1.0e+00), Tii, ldt, T0i, ldt );

	}

	return 0;

}


