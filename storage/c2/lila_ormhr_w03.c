#include "lila.h"

int lila_ormhr_w03( int m, int n, int i, int j, int l, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork, int *S ){

	double *Q0j, *Qij, *Qjj, *Tij, *Tjj, *Ajj, *Aii, *Aij, *Aji;
	int k, i1, info;

	Q0j = Q + j*ldq;
	Qij = Q + i + j*ldq;
	Qjj = Q + j + j*ldq;

	Tij = T + (i % mt) + j*ldt;
	Tjj = T + (j % mt) + j*ldt;

	Aij = A + i + j*lda;
	Aii = A + i + i*lda;
	Ajj = A + j + j*lda;
	Aji = A + j + i*lda;

	double *Ali, *All, *Ajl, *Qlj, *Tlj;
	Tlj = T + (l % mt) + j*ldt;
	Ajl = A + j + l*lda;
	Ali = A + l + i*lda;
	All = A + l + l*lda;
	Qlj = Q + l + j*ldq;

	for( k = 0; k < j; k++){
		if( S[k] == -1 ){
			for( i1 = 0; i1 < n; i1++) Aij[ k + i1*lda ] = - Aij[ k + i1*lda ];
		}
	}

	double *work2;
	int lwork2;

	lwork2 = j-l;
	work2 = (double *) malloc(lwork2 * n * sizeof(double));

//	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', j-l, n, Qlj, ldq, Tlj, ldt ); 
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', j-l, n, Qlj, ldq, work2, lwork2 ); 

	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', l-i, n, Qij, ldq, work, lwork ); 

	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', n, n, Qjj, ldq, Tjj, ldt ); 
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m-j-i-n, n, Qjj+n, ldq, Ajj+n, lda );

	if( l > i ) cblas_dtrsm( CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, l-i, n, 1.0e+00, Aii, lda, work, lwork );

//	if( l > i ) cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, j-l, n, l-i, (-1.0e+00), Ali, lda, work, lwork, (+1.0e+00), Tlj, ldt );
	if( l > i ) cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, j-l, n, l-i, (-1.0e+00), Ali, lda, work, lwork, (+1.0e+00), work2, lwork2 );

//	cblas_dtrsm( CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, j-l, n, 1.0e+00, All, lda, Tlj, ldt );
	if( lwork2 != 0 ) cblas_dtrsm( CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, j-l, n, 1.0e+00, All, lda, work2, lwork2 );

	if( l > i ) cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, l-i, (-1.0e+00), Aji, lda, work, lwork, (+1.0e+00), Tjj, ldt );
 
//	cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, j-l, (-1.0e+00), Ajl, lda, Tlj, ldt, (+1.0e+00), Tjj, ldt ); 
	cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, j-l, (-1.0e+00), Ajl, lda, work2, lwork2, (+1.0e+00), Tjj, ldt ); 

	if( l > i ) cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, m-j-i-n, n, l-i, (-1.0e+00), Aji+n, lda, work, lwork, (+1.0e+00), Ajj+n, lda ); 
 
//	cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, m-j-i-n, n, j-l, (-1.0e+00), Ajl+n, lda, Tlj, ldt, (+1.0e+00), Ajj+n, lda ); 
	cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, m-j-i-n, n, j-l, (-1.0e+00), Ajl+n, lda, work2, lwork2, (+1.0e+00), Ajj+n, lda ); 

	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', l-i, n, work, lwork, Tij, ldt ); 

	free( work2 );
	
	return 0;

}

