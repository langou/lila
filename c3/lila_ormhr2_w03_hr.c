#include "lila.h"

int lila_ormhr2_w03_hr( int m, int n, int i, int j, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork, int *S ){

	double *Qij, *Qjj, *Tjj, *Tkk, *Ajj, *Aii, *Aij, *Aji, *A0j, *work2;
	int k, info, i1, it, lwork2;

	it = (j % mt); 
	lwork2 = j-i;
//	work2 = (double *) malloc(lwork2 * n * sizeof(double));
	work2 = work;

	Qij = Q + i  + j*ldq;
	Qjj = Q + j  + j*ldq;

	Tjj = T + it + j*ldt;

	A0j = A +      j*lda;
	Aii = A + i  + i*lda;
	Ajj = A + j  + j*lda;
	Aij = A + i  + j*lda;
	Aji = A + j  + i*lda;

	for( k = i; k < j; k++){
		if( S[ k ] == -1 ){
			for( i1 = 0; i1 < n; i1++) A0j[ k + i1*lda ] = - A0j[ k + i1*lda ];
		}
	}

//	printf("lwork = %d\n\n",lwork);

	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A',   j-i, n,   Qij, ldq, work2, lwork2 ); 
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A',     n, n,   Qjj, ldq,   Tjj,    ldt ); 
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m-j-n, n, Qjj+n, ldq, Ajj+n,    lda ); 

	if( j-i != 0) cblas_dtrsm( CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, j-i, n, (+1.0e+00), Aii, lda, work2, lwork2 );
	if( j-i != 0) cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans,     n, n, j-i, (-1.0e+00),   Aji, lda, work2, lwork2, (+1.0e+00),   Tjj, ldt ); 
	if( j-i != 0) cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, m-j-n, n, j-i, (-1.0e+00), Aji+n, lda, work2, lwork2, (+1.0e+00), Ajj+n, lda ); 

	for( k = 0; k < n; k++ ){

		Tkk = Tjj + k + k*ldt;
		if ( fabs( 1.0e+00 - (*Tkk) ) < fabs( 1.0e+00 + (*Tkk) ) ){
			for( i1 = 0; i1 < n; i1++) Tjj[ i1 + k*ldt ] = - Tjj[ i1 + k*ldt];
 	 		S[ k+j ] = - 1;
		} else {
			S[ k+j ] =   1;		
		}

		(*Tkk) = (*Tkk) - 1.0e+00;
		for( i1 = k+1; i1 < n; i1++ ) Tjj[ i1 + k*ldt ]  = Tjj[ i1 + k*ldt ] / (*Tkk);
		cblas_dger( CblasColMajor, n-k-1, n-k-1, (-1.0e+00), Tkk+1, 1, Tkk+ldt, ldt, Tkk+1+ldt, ldt );

	}

	for( k = 0; k < n; k++ ){
		if ( S[ k+j ] == -1 ){
			for( i1 = k; i1 < n; i1++) Ajj[ k + i1*lda ] = - Ajj[ k + i1*lda ];
		}
	}

	for( k = 0; k < n; k++ ){
		if ( S[ k+j ] == -1 ){
			for( i1 = n; i1 < m-j; i1++) Ajj[ i1 + k*lda ] = - Ajj[ i1 + k*lda ]; 
		}
	}

	for( k = 0; k < n; k++ ){
		if ( S[ k+j ] == -1 ){
			for( i1 = 0; i1 < m-i; i1++) Qij[ i1 + k*ldq ] = - Qij[ i1 + k*ldq ]; 
		}
	}

	for( k = 0; k < n; k++ ){
		for( i1 = 0; i1 < k; i1++) Ajj[ k + i1*lda ] = Tjj[ k + i1*ldt ];
	}

	for( k = 0; k < n; k++ ){
		for( i1 = 0; i1 < k; i1++) Tjj[ k + i1*ldt ] = (+0.0e+00);
	}

	cblas_dtrsm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, m-j-n, n, (+1.0e+00), Tjj, ldt, Ajj+n, lda ); 
	cblas_dtrsm( CblasColMajor, CblasRight, CblasLower,   CblasTrans,    CblasUnit,     n, n, (-1.0e+00), Ajj, lda,   Tjj, ldt );

//	free( work2 );

	return 0;

}

