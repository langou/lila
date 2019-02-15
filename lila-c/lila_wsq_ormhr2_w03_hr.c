#include "lila.h"

int lila_wsq_ormhr2_w03_hr( int m, int n, int i, int j, int l, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork, double *S ){

//	lwork2 = j-i;
//	work2  = (double *) malloc(lwork2 * n * sizeof(double));
//	work2 = work + n+i;

//	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A',   j-i, n,   Qij, ldq, work2, lwork2 ); 

//	if( j-i != 0) cblas_dtrsm( CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, j-i, n, (+1.0e+00), Aii, lda, work2, lwork2 );
//	if( j-i != 0) cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans,     n, n, j-i, (-1.0e+00),   Aji, lda, work2, lwork2, (+1.0e+00),   Tjj, ldt ); 
//	if( j-i != 0) cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, m-j-n, n, j-i, (-1.0e+00), Aji+n, lda, work2, lwork2, (+1.0e+00), Ajj+n, lda ); 

	return (n+i)*(j-i); // These are the dimension that I set-up the array. But do I need that or just j-i?
//	return (j-i); // These are the dimension that I set-up the array. But do I need that or just j-i?

}

