#include "V2T.h"

int UinvLTinU( int n, double *L, int ldl, double *U, int ldu ){

	int n1, n2, i, j; 
	double *U11, *U12, *U22, *L11, *L21, *L22;

	if ( n <= 1 ){

		(*U) = 1.0e+00 / (*U);

	} else { 

		n1 = n/2; n2 = n-n1;

		L11 = L;
		U11 = U;

		U12 = U + n1*ldu;
		L21 = L + n1; 

		U22 = U + n1 + n1*ldu;  
		L22 = L + n1 + n1*ldl; 

		UinvLTinU( n2, L22, ldl, U22, ldu );

		cblas_dtrmm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, n1, n2, (-1.0e+00), U22, ldu, U12, ldu );
		for(i=0;i<n1;i++){for(j=0;j<n2;j++){ U12[i+j*ldu] += L21[j+i*ldl]; }}
		cblas_dtrsm( CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, n1, n2, (+1.0e+00), U11, ldu, U12, ldu );			

		UinvLTinU( n1, L11, ldl, U11, ldu );

	}

	return 0;

}
