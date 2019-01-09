#include "lila.h"

int lila_dorgh2_3( int m, int n, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork, int* S ){

	double *Tkk;
	int k, i1;

	for( k = 0; k < n; k++ ){
		Tkk = T + k + k*ldt;

		if ( fabs( 1.0e+00 - (*Tkk) ) < fabs( 1.0e+00 + (*Tkk) ) ){
			for( i1 = 0; i1 < n; i1++) T[ i1 + k*ldt ] = - T[ i1 + k*ldt];
 	 		S[ k ] = - 1;
		} else {
			S[ k ] = 1;		
		}
		(*Tkk) = (*Tkk) - 1.0e+00;
		for( i1 = k+1; i1 < n; i1++ ) T[ i1 + k*ldt ]  = T[ i1 + k*ldt ] / (*Tkk);
		cblas_dger( CblasColMajor, n-k-1, n-k-1, (-1.0e+00), Tkk+1, 1, Tkk+ldt, ldt, Tkk+1+ldt, ldt );

	}

	for( k = 0; k < n; k++ ){
		if ( S[ k ] == -1 ){
			for( i1 = n; i1 < m; i1++) A[ i1 + k*lda ] = - A[ i1 + k*lda ]; 
		}
	}

	cblas_dtrsm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, m-n, n, 1.0e+00, T, ldt, A+n, lda ); 

	for( k = 0; k < n; k++ ){
		if ( S[ k ] == -1 ){
			for( i1 = k; i1 < n; i1++) A[ k + i1*lda ] = - A[ k + i1*lda ];
		}
	}

	for( k = 0; k < n; k++ ){
		if ( S[ k ] == -1 ){
			for( i1 = 0; i1 < m; i1++) Q[ i1 + k*ldq ] = - Q[ i1 + k*ldq ]; 
		}
	}

	for( k = 0; k < n; k++ ){
		for( i1 = 0; i1 < k; i1++) A[ k + i1*lda ] = T[ k + i1*ldt ];
	}

	for( k = 0; k < n; k++ ){
		for( i1 = 0; i1 < k; i1++) T[ k + i1*ldt ] = 0.0e+00;
	}

	cblas_dtrsm( CblasColMajor, CblasRight, CblasLower, CblasTrans, CblasUnit, n, n, (-1.0e+00), A, lda, T, ldt );

	return 0;

}
