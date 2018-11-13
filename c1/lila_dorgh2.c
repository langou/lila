#include "lila.h"

int lila_dorgh2( int m, int n, int j, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork, int* S ){

	double *Tjj, *Ajj, *T0j;
	double *Tkk;
	int *Sj;
	int k, i1;

	Tjj = T + j + j*ldt;
	T0j = T + j*ldt;
	Ajj = A + j + j*lda;
	Sj = S + j;

//	printf("\n j = %d, n = %d, m = %d,\n",j,n,m);

	for( k = 0; k < n; k++ ){

		Tkk = Tjj + k + k*ldt;

		if ( fabs( 1 - Tjj[ k + k*ldt ] ) < fabs( 1 + Tjj[ k + k*ldt ] ) ){

			for( i1 = 0; i1 < n; i1++) Tjj[ i1 + k*ldt ] = - Tjj[ i1 + k*ldt];
 	 		Sj[ k ] = - 1;		

		} else {
			
			Sj[ k ] = - 1;		
	
		}
		
		Tjj[ k + k*ldt ] -= 1;
		for( i1 = k; i1 < n; i1++ ) Tjj[ i1 + k*ldt ]  = Tjj[ i1 + k*ldt ] / Tjj[ k + k*ldt ];
		cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, n-k, n-k, 1, (-1.0e+00), Tkk, ldt, Tkk, ldt, (+1.0e+00), Tjj, ldt ); 

	}		

	for( k = 0; k < n; k++ ){
		
		if ( Sj[ k ] == -1 ){

			for( i1 = n; i1 < m; i1++) Ajj[ i1 + k*lda ] = - Ajj[ i1 + k*lda ];

		}
	}

	cblas_dtrsm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, m-j-n, n, 1.0e+00, Tjj, ldt, Ajj + n, lda );

	for( k = 0; k < n; k++ ){
		
		if ( Sj[ k ] == -1 ){

			for( i1 = 0; i1 < j; i1++) T0j[ i1 + k*ldt ] = - T0j[ i1 + k*ldt ];

		}
	}

	for( k = 0; k < n; k++ ){
		
		if ( Sj[ k ] == -1 ){

			for( i1 = k; i1 < n; i1++) Ajj[ k + i1*lda ] = - Ajj[ k + i1*lda ];

		}
	}

	for( k = 0; k < n; k++ ){
		
		if ( Sj[ k ] == -1 ){

			for( i1 = 0; i1 < m; i1++) Q[ i1 + (j+k)*ldq ] = - Q[ i1 + (j+k)*ldq ];

		}
	}

	for( k = 0; k < n; k++ ){

		for( i1 = 0; i1 < k; i1++) Ajj[ k + i1*lda ] = Tjj[ k + i1*ldt ];

	}

	cblas_dtrsm( CblasColMajor, CblasRight, CblasLower, CblasTrans, CblasUnit, n, n, (-1.0e+00), Ajj, lda, Tjj, ldt );



	return 0;














// This potion below is what we had before I started the day
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	this subroutine needs a workspace of size n-by-n
//	it is possible to reduce this workspace to nothing by changing the algorithm
//	but we think it is fine like this
/*
	int info, i1, j1;
	int iiii;


//	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'U', n, n, Q, ldq, T, ldt );
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'U', n, n, T, ldt, work, n );
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'L', n-1, n, A+1, lda, work+1, n );

	for( i1 = 0; i1 < n ; i1++){

//		if ( fabs(1 - work[ i1 + i1*n ]) < fabs(1 + work[ i1 + i1*n ]) ){

	
			//S[ i1 ] = -1;
//
		//} else {

			S[ i1 ] = 1;
			for( j1 = 0; j1 < n ; j1++) work[ j1 + i1*n ] = - work[ j1 + i1*n ];			
			for( j1 = n; j1 < m ; j1++) A[ j1 + i1*lda ] = - A[ j1 + i1*lda ];			

		//}

		work[ i1 + i1*n ] += 1.0e+00; 

		for( j1 = i1+1; j1 < n; j1++ ) work[ j1 + i1*n ] = work[ j1 + i1*n ] / work[ i1 + i1*n ];
		for( j1 = n; j1 < m; j1++ ) A[ j1 + i1*lda ] = A[ j1 + i1*lda ] / work[ i1 + i1*n ];

		cblas_dger( CblasColMajor, n-i1-1, n-i1-1, (-1.0e+00), work + (i1+1) + i1*n, 1, work + i1 + (i1+1)*n, n, work + (i1+1) + (i1+1)*n, n);
		cblas_dger( CblasColMajor, m-n, n-i1-1, (-1.0e+00), A + n + i1*lda, 1, work + i1 + (i1+1)*n, n, A + n + (i1+1)*lda, lda);

	}

	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'U', n, n, work, n, T, ldt ); 
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'L', n-1, n, work+1, n, A+1, lda );

	cblas_dtrsm( CblasColMajor, CblasRight, CblasLower, CblasTrans, CblasUnit, n, n, 1.0e+00, work, n, T, ldt );

	return 0;
*/

}







