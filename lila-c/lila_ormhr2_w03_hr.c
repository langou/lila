#include "lila.h"

int lila_ormhr2_w03_hr( int m, int n, int i, int j, int l, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork, int *S ){

	double *Q0j, *Qij, *Qjj, *Tij, *Tjj, *Ajj, *Aii, *Aij, *Aji, *Ai0, *Aj0;
	int k, info, i1;

	Q0j = Q + j*ldq;
	Qij = Q + i + j*ldq;
	Qjj = Q + j + j*ldq;

	Tij = T + (i % mt) + j*ldt;
	Tjj = T + (j % mt) + j*ldt;

	Ai0 = A + i;
	Aj0 = A + j*lda;
	Aij = A + i + j*lda;
	Aii = A + i + i*lda;
	Ajj = A + j + j*lda;
	Aji = A + j + i*lda;

	for( k = 0; k < j; k++){
		if( S[k] == -1 ){
			for( i1 = 0; i1 < n; i1++) Aij[ k + i1*lda ] = - Aij[ k + i1*lda ];
		}
	}

	double *work2;
	int lwork2;

	lwork2 = j-i;
	work2 = (double *) malloc(lwork2 * n * sizeof(double));

	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', j-i, n, Qij, ldq, work2, lwork2 ); 

	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', n, n, Qij, ldq, Tij, ldt ); 

	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m-j-i-n, n, Qjj+n, ldq, Ajj+n, lda ); 

//	if( i != 0 ) cblas_dtrsm( CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, i, n, 1.0e+00, A, lda, work2, lwork2 );

//	if( i != 0 ) cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, j, n, i, (-1.0e+00), Aii, lda, work2, lwork2, (+1.0e+00), Tij, ldt );

	cblas_dtrsm( CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, j-i, n, 1.0e+00, Aii, lda, Tij, ldt );

//	if( i != 0 ) cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, i, (-1.0e+00), Aji, lda, work2, lwork2, (+1.0e+00), Tjj, ldt );

	cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, j, (-1.0e+00), Aji, lda, Tij, ldt, (+1.0e+00), Tjj, ldt ); 

//	if( i != 0 ) cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, m-j-i-n, n, i, (-1.0e+00), Aji+n, lda, work2, lwork2, (+1.0e+00), Ajj+n, lda ); 

	cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, m-j-i-n, n, j-i, (-1.0e+00), Aji+n, lda, Tij, ldt, (+1.0e+00), Ajj+n, lda ); 

// we   info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', l-i, n, work, lwork, Tij, ldt ); 
// do	if( i != 0 ) info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', i, n, work, lwork, Tij, ldt ); 
// not
// have in matlab script 

	printf("\n\n");

	free( work2 );
	
	return 0;

}

