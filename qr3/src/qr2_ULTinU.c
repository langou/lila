#include "qr2.h"

int qr2_ULTinU( int n, double *L, int ldl, double *U, int ldu ){

	int n1, n2; 
	double *U11, *U12, *U22, *L11, *L21, *L22;

	if ( n <= 1 ){

		//Nothing needs to be done

	} else { 

		n1 = n/2; n2 = n-n1;

		L11 = L;
		U11 = U;

		U12 = U + n1*ldu;
		L21 = L + n1; 

		U22 = U + n1 + n1*ldu;  
		L22 = L + n1 + n1*ldl; 

//		U12 = U12*L22^T
		cblas_dtrmm( CblasColMajor, CblasRight, CblasLower, CblasTrans, CblasUnit, n1, n2, (+1.0e+00), L22, ldl, U12, ldu );

//		U12 = U12 + U11*L21^T
		qr2_ApUBTinA ( n1, n2, U12, ldu, U11, ldu, L21, ldl );

		qr2_ULTinU( n1, L11, ldl, U11, ldu );
		qr2_ULTinU( n2, L22, ldl, U22, ldu );
		
	}

	return 0;

}
