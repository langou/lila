#include "lila.h"

int lila_dge_qr_wq_levelx_w02( int n_lvl, int i_lvl, int *nb_lvl, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork ){

	int vb, info; 
	int j, k, nb;

	nb = nb_lvl[ i_lvl ];

	if ( nb > n ) vb = n; else vb = nb; 

	k = 0;
	j = i;

	if( i_lvl == n_lvl-1 ){
	printf("vb (level %d) = %d\n", i_lvl, vb );
	info = lila_dge_qr_wq_w02( m, vb, j, mt, A, lda, T, ldt, Q, ldq, work, lwork );
	} else {
	printf("vb (level %d) = %d\n", i_lvl, vb );
	info = lila_dge_qr_wq_levelx_w02( n_lvl, i_lvl+1, nb_lvl, m, vb, j, mt, A, lda, T, ldt, Q, ldq, work, lwork );
	}

	k += vb;
	j += vb;
	
	while ( j < i+n ) {

	if ( j+nb > i+n ) vb = i+n-j; else vb = nb; 

	info = lila_dge_qr_ormqrf_w02( m, vb, k, i, j, mt, A, lda, T, ldt, work, lwork );

	if( i_lvl == n_lvl-1 ){
	printf("vb (level %d) = %d\n", i_lvl, vb );
	info = lila_dge_qr_wq_w02( m, vb, j, mt, A, lda, T, ldt, Q, ldq, work, lwork );
	} else {
	printf("vb (level %d) = %d\n", i_lvl, vb );
	info = lila_dge_qr_wq_levelx_w02( n_lvl, i_lvl+1, nb_lvl, m, vb, j, mt, A, lda, T, ldt, Q, ldq, work, lwork );
	}

	int ii, jj;
	double *Aii, *Aij, *Aji, *Ajj, *Tii, *Tij, *Tjj;

	Aii = A + i + i*lda;
	Aij = A + i + j*lda;
	Aji = A + j + i*lda;
	Ajj = A + j + j*lda;

	Tii = T + i + i*ldt;
	Tij = T + i + j*ldt;
	Tjj = T + j + j*ldt;

	for( jj = 0; jj < vb; jj++ ){
	for( ii = 0; ii < k; ii++ ){
	Tij[ ii + jj * ldt ] = Aji[  jj + ii * lda  ];
	}}

	cblas_dtrmm ( CblasColMajor, CblasRight, CblasLower, CblasNoTrans, CblasUnit, k, vb, (+1.0e+00), Ajj, lda, Tij, ldt );
	cblas_dgemm ( CblasColMajor, CblasTrans, CblasNoTrans, k, vb, m-j-vb, (+1.0e+00), Aji+vb, lda, Ajj+vb, lda, (+1.0e+00), Tij, ldt );
	cblas_dtrmm ( CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, k, vb, (-1.0e+00), Tii, ldt, Tij, ldt );
	cblas_dtrmm ( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, k, vb, (+1.0e+00), Tjj, ldt, Tij, ldt );

	info = lila_dge_qr_ormqrbz_w02( m, vb, k, i, j, mt, A, lda, Q, ldq, T, ldt, work, lwork );
 
	k += vb;
	j += vb;

	}

	return 0;

}
