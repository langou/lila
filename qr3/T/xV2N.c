#include "lila.h"

int xV2N( int n, double *T, int ldt ){

	//double *Tii, *taui;
	//double *T0i;
	//int info, n1, n2, ii, jj; 

	double *V; int i,j;
	V  = (double *) malloc( n * n * sizeof(double));
	cblas_dsyrk( CblasColMajor, CblasUpper, CblasNoTrans, n, n, (+1.0e+00), T, ldt, (+0.0e+00), V, n );

	for(i=0;i<n;i++){for(j=i;j<n;j++){ T[i+j*ldt] = V[i+j*n];}}

	free( V );

/*

	if ( n <= 1 ) {

//		(*T) = (*tau);
		(*T) = (+2.0e+00) / (*T);

	} else {

		n1 = n / 2;
		n2 = n - n1;

		Tii = T + n1 + n1*ldt;
		taui = tau + n1;

		T0i = T + n1*ldt;

		info = xV2N( n1, tau, T, ldt );
		info = xV2N( n2, taui, Tii, ldt );

		cblas_dtrmm( CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit,  n1, n2, (-1.0e+00), T, ldt, T0i, ldt );
		cblas_dtrmm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, n1, n2, (+1.0e+00), Tii, ldt, T0i, ldt );

	}

	return 0;
*/

}


