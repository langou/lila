#include "lila.h"

int lila_ormhr_w03( int m, int n, int i, int j, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, int *S ){

	double *work, *Q0j, *Qij, *Qjj, *Tij, *Tjj, *Ajj, *Aii, *Aij, *Aji;
	int k, l, i1, info, ml, ldw, vb;

//	vb = mt - (j % mt ); if ( vb > n ) vb = n;
	ml = m-i;

	l = 5;
	if ( j > i + l  ) l = i + l; else l = i;

	ldw = l-i;
	work = (double *) malloc(ldw * n * sizeof(double));

	Q0j = Q + j*ldq;
	Qij = Q + i + j*ldq;
	Qjj = Q + j + j*ldq;

	Tij = T + (i % mt) + j*ldt;
	//Tij = T + i + j*ldt;
	Tjj = T + (j % mt) + j*ldt;
	//Tjj = T + j + j*ldt;

	Aij = A + i + j*lda;
	Aii = A + i + i*lda;
	Ajj = A + j + j*lda;
	Aji = A + j + i*lda;

	double *Tlj, *Qlj;
	Tlj = T + (l % mt) + j*ldt;
	//Tlj = T + l + j*ldt;
	Qlj = Q + l + j*ldq;

	double *Ali, *All, *Ajl, *Til;
	Ajl = A + j + l*lda;
	Ali = A + l + i*lda;
	All = A + l + l*lda;
	Til = T + (i % mt) + l*ldt;
	//Til = T + i + l*ldt;

	for( k = 0; k < j; k++){
	
		if( S[k] == -1 ){

			for( i1 = 0; i1 < n; i1++) Aij[ k + i1*lda ] = - Aij[ k + i1*lda ];

		}

	}

//	int not_done, jj;
//	jj = 1;
//	not_done = 1;

//	while( not_done == 1 ){

	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', j-l, n, Qlj, ldq, Tlj, ldt ); 
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', l-i, n, Qij, ldq, work, ldw ); 

	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', n, n, Qjj, ldq, Tjj, ldt ); 
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml-j-n, n, Qjj+n, ldq, Ajj+n, lda ); 

	if( l > i ) cblas_dtrsm( CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, l-i, n, 1.0e+00, Aii, lda, work, ldw );
	if( l > i ) cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, j-l, n, l-i, (-1.0e+00), Ali, lda, work, ldw, (+1.0e+00), Tlj, ldt ); 
	cblas_dtrsm( CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, j-l, n, 1.0e+00, All, lda, Tlj, ldt );

	if( l > i ) cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, l-i, (-1.0e+00), Aji, lda, work, ldw, (+1.0e+00), Tjj, ldt ); 
	cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, j-l, (-1.0e+00), Ajl, lda, Tlj, ldt, (+1.0e+00), Tjj, ldt ); 

	if( l > i ) cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, ml-j-n, n, l-i, (-1.0e+00), Aji+n, lda, work, ldw, (+1.0e+00), Ajj+n, lda ); 
	cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, ml-j-n, n, j-l, (-1.0e+00), Ajl+n, lda, Tlj, ldt, (+1.0e+00), Ajj+n, lda ); 

	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', l-i, n, work, ldw, Tij, ldt ); 

//		if ( jj + vb - 1 == n ) {
//			not_done = 0; 
//		} else if ( jj + vb - 1 > n ) {
//			printf("defensive programming, we should never have been there, abort\n"); return 0;
//		} else {
//			if(jj == 1) Tjj = Tjj - (j%mt) + vb*ldt; 
//			ml -= vb;
//			jj += vb;
//			if(jj != 1+vb) Tjj = Tjj + vb*ldt; 
//			if ( ( jj + mt - 1 ) <= n ) vb = mt; else vb = n - jj + 1;
//		}
//	}	
	free( work );

	return 0;

}

