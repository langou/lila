#include "lila.h"

int lila_dorgh2( int m, int n, int i, int j, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork, int* S ){

	double *Tjj, *Ajj, *T0j;
	double *Tkk, *Q0j;
	int *Sj;
	int k, i1, ml;

	ml = m;

	Tjj = T + j + j*ldt;
	T0j = T + j*ldt;
	Ajj = A + j + j*lda;
	Q0j = Q + j*ldq;
	Sj = S + j;

	printf("\n j = %d, n = %d, m = %d,\n",j,n,m);

//	LU Factorization
	for( k = 0; k < n; k++ ){
		Tkk = Tjj + k + k*ldt;

		if ( fabs( 1.0e+00 - (*Tkk) ) < fabs( 1.0e+00 + (*Tkk) ) ){
			for( i1 = 0; i1 < n; i1++) Tjj[ i1 + k*ldt ] = - Tjj[ i1 + k*ldt];
 	 		Sj[ k ] = - 1;		
			//Sj[ k ] = 1;		
		} else {
			Sj[ k ] = 1;		
		}
		(*Tkk) = (*Tkk) - 1.0e+00;
		for( i1 = k+1; i1 < n; i1++ ) Tjj[ i1 + k*ldt ]  = Tjj[ i1 + k*ldt ] / (*Tkk);
		cblas_dger( CblasColMajor, n-k-1, n-k-1, (-1.0e+00), Tkk+1, 1, Tkk+ldt, ldt, Tkk+1+ldt, ldt );

	}

	for( k = 0; k < n; k++ ){
		if ( Sj[ k ] == -1 ){
			for( i1 = n; i1 < m-j; i1++) Ajj[ i1 + k*lda ] = - Ajj[ i1 + k*lda ];
		}
	}

	cblas_dtrsm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, m-j-n, n, 1.0e+00, Tjj, ldt, Ajj + n, lda );

	for( k = 0; k < n; k++ ){
		if ( Sj[ k ] == -1 ){
			for( i1 = 0; i1 < j-1; i1++) T0j[ i1 + k*ldt ] = - T0j[ i1 + k*ldt ];
		}
	}

	for( k = 0; k < n; k++ ){
		if ( Sj[ k ] == -1 ){
			for( i1 = k; i1 < n; i1++) Ajj[ k + i1*lda ] = - Ajj[ k + i1*lda ];
		}
	}

	for( k = 0; k < n; k++ ){
		if ( Sj[ k ] == -1 ){
			for( i1 = 0; i1 < m-i; i1++) Q0j[ i1 + k*ldq ] = - Q0j[ i1 + k*ldq ];
		}
	}

	for( k = 0; k < n; k++ ){
		for( i1 = 0; i1 < k; i1++) Ajj[ k + i1*lda ] = Tjj[ k + i1*ldt ];
	}

	for( k = 0; k < n; k++ ){
		for( i1 = 0; i1 < k; i1++) Tjj[ k + i1*ldt ] = 0.0e+00;
	}

	cblas_dtrsm( CblasColMajor, CblasRight, CblasLower, CblasTrans, CblasUnit, n, n, (-1.0e+00), Ajj, lda, Tjj, ldt );

	return 0;

}
