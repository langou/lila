#include "lila.h"

int lila_dorgh2( int m, int n, int j, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork, int* S ){


	double *Tjj, *Ajj, *T0j;
	double *Tkj, *Tjk;
	int *Sj;
	int k, i1;

	Tjj = T + j + j*ldt;
	T0j = T + j*ldt;
	Ajj = A + j + j*lda;
	Sj = S + j;

	for( k = 0; k < n; k++ ){

		Tkj = Tjj + k;
		Tjk = Tjj + k*ldt;

		if ( fabs( 1 - Tjj[ k + k*ldt ] ) < fabs( 1 + Tjj[ k + k*ldt ] ) ){

			//for( i1 = 0; i1 < n; i1++) Tjj[ i1 + k*ldt ] = - Tjj[ i1 + k*ldt];
 	 		Sj[ k ] = - 1;		

		} else {
			
			Sj[ k ] = - 1;		
	
		}
		
		//Tjj[ k + k*ldt ] -= 1;
		//for( i1 = k; i1 < n; i1++ ) Tjj[ i1 + k*ldt ]  = Tjj[ i1 + k*ldt ] / Tjj[ k + k*ldt ];
		//cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, n-k, n-k, 1, (-1.0e+00), Tkj, ldt, Tjk, ldt, (+1.0e+00), Tjj, ldt ); 

	}		

	for( k = 0; k < n; k++ ){
		
		if ( Sj[ k ] == -1 ){

			//for( i1 = n; i1 < m; i1++) Ajj[ i1 + k*lda ] = - Ajj[ i1 + k*lda ];

		}
	}

	//cblas_dtrsm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, m-j-n, n, 1.0e+00, Tjj, ldt, Ajj + n, lda );

	for( k = 0; k < n; k++ ){
		
		if ( Sj[ k ] == -1 ){

			//for( i1 = 0; i1 < j; i1++) T0j[ i1 + k*ldt ] = - T0j[ i1 + k*ldt ];

		}
	}

	for( k = 0; k < n; k++ ){
		
		if ( Sj[ k ] == -1 ){

			//for( i1 = k; i1 < n; i1++) Ajj[ k + i1*lda ] = - Ajj[ k + i1*lda ];

		}
	}

	for( k = 0; k < n; k++ ){
		
		if ( Sj[ k ] == -1 ){

			//for( i1 = 0; i1 < m; i1++) Q[ i1 + (j+k)*ldq ] = - Q[ i1 + (j+k)*ldq ];

		}
	}

	for( k = 0; k < n; k++ ){

		//for( i1 = 0; i1 < k; i1++) Ajj[ k + i1*lda ] = Tjj[ k + i1*ldt ];

	}

	//cblas_dtrsm( CblasColMajor, CblasRight, CblasLower, CblasTrans, CblasUnit, n, n, (-1.0e+00), Ajj, lda, Tjj, ldt );



	return 0;
}







