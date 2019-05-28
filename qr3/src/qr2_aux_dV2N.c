#include "qr2.h"

int qr2_aux_dV2N( int n, double *T, int ldt ){

	int n1, n2;
	double *T11, *T12, *T22;

	if ( n <= 1 ) {

		(*T) = 1.0e+00;

	} else {

		n1 = n / 2;
		n2 = n - n1;

		T11 = T;
		T12 = T+n1*ldt;
		T22 = T+n1+n1*ldt;

		qr2_aux_dV2N( n1, T11, ldt );
		cblas_dsyrk( CblasColMajor, CblasUpper, CblasNoTrans, n1, n2, (+1.0e+00), T12, ldt, (+1.0e+00), T11, ldt );
		cblas_dtrmm( CblasColMajor, CblasRight, CblasUpper, CblasTrans, CblasUnit,  n1, n2, (+1.0e+00), T22, ldt, T12, ldt );
		qr2_aux_dV2N( n2, T22, ldt );

	}

	return 0;

}


