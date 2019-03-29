#include "lila.h"

int ULTinU( int n, double *L, int ldl, double *U, int ldu ){

	int info, n1, n2; 
	double *work, *U11, *U12, *U22, *L11, *L12, *L22;

	if ( n <= 1 ){

		//Nothing needs to be done

	} else { 

		n1 = n/2; n2 = n-n1;

		L11 = L;
		U11 = U;

		U12 = U + n1*ldu;
		L12 = L + n1; 

		U22 = U + n1 + n1*ldu;  
		L22 = L + n1 + n1*ldl; 

		cblas_dtrmm( CblasColMajor, CblasRight, CblasLower, CblasTrans, CblasUnit, n1, n2, (+1.0e+00), L22, ldl, U12, ldu );
		cblas_dgemm( CblasColMajor, CblasNoTrans, CblasTrans, n1, n2, n1, (+1.0e+00), U11, ldu, L12, ldl, (+1.0e+00), U12, ldu );

		ULTinU( n1, L11, ldl, U11, ldu );
		ULTinU( n2, L22, ldl, U22, ldu );
		
	}

	return 0;

}
