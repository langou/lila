#include "lila.h"

int lila_ormhr_w0b( int m, int n, int i, int j, double *A, int lda, double *T, int ldt, double *Q, int ldq, int *S ){

	double *work, *Q0j, *Qij, *Qjj, *Tij, *Tjj, *Ajj, *Ai0, *Aii, *Aj0, *A0j, *Aji;
	int k, i1, info;

	//work = (double *) malloc((i+n) * n * sizeof(double));
	work = (double *) malloc(m * n * sizeof(double));

	Q0j = Q + j*ldq;
	Qij = Q + i + j*ldq;
	Qjj = Q + j + j*ldq;

	Tij = T + i + j*ldt;
	Tjj = T + j + j*ldt;

	Ai0 = A + i;
	Aj0 = A + j;
	A0j = A + j*lda;
	Aii = A + i + i*lda;
	Ajj = A + j + j*lda;
	Aji = A + j + i*lda;


	for( k = 0; k < j-1; k++){
	
		if( S[k] == -1 ){

			for( i1 = 0; i1 < n; i1++) A0j[ k + i1*lda ] = - A0j[ k + i1*lda ];

		}

	}

	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', i, n, Q0j, ldq, work, n ); 
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', j-i, n, Qij, ldq, Tij, ldt ); 
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', n, n, Qjj, ldq, Tjj, ldt ); 
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m-j-n+1, n, Qjj+n, ldq, Ajj+n, lda ); 

	cblas_dtrsm( CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, i, n, 1.0e+00, A, lda, work, n );
	cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, j-i, n, i, (-1.0e+00), Ai0, lda, work, n, (+1.0e+00), Tij, ldt ); 
	cblas_dtrsm( CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, j-i, n, 1.0e+00, Aii, lda, Tij, ldt );
	cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, i, (-1.0e+00), Aj0, lda, work, n, (+1.0e+00), Tjj, ldt ); 
	cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, j-i, (-1.0e+00), Aji, lda, Tij, ldt, (+1.0e+00), Tjj, ldt ); 
	cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, m-j-n, n, i, (-1.0e+00), Aj0+n, lda, work, n, (+1.0e+00), Ajj+n, lda ); 
	cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, m-j-n, n, j-i, (-1.0e+00), Aji+n, lda, Tij, ldt, (+1.0e+00), Ajj+n, lda ); 

	free( work );

	return 0;

}
