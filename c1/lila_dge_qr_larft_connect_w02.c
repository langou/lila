#include "lila.h"

//int lila_dge_qr_larft_connect_w02( int m, int n, int k, int i, int j, int mt, double *A, int lda, double *T, int ldt ){
int lila_dge_qr_larft_connect_w02( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt ){


	double *Aii, *Ai0, *Tii, *T0i;
	int ii, jj;

	Aii = A + i + i*lda;
	Tii = T + i + i*ldt;
	T0i = T + i*ldt;
	Ai0 = A + i;

	for( jj = 0; jj < n; jj++ ){
	for( ii = 0; ii < i; ii++ ){
		T0i[ ii + jj * ldt ] = Ai0[  jj + ii * lda  ];
	}}

	cblas_dtrmm ( CblasColMajor, CblasRight, CblasLower, CblasNoTrans, CblasUnit, i, n, (+1.0e+00), Aii, lda, T0i, ldt );
	cblas_dgemm ( CblasColMajor, CblasTrans, CblasNoTrans, i, n, m-n-i, (+1.0e+00), Ai0+n, lda, Aii+n, lda, (+1.0e+00), T0i, ldt );
	cblas_dtrmm ( CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, i, n, (-1.0e+00), T, ldt, T0i, ldt );
	cblas_dtrmm ( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, i, n, (+1.0e+00), Tii, ldt, T0i, ldt );


	return 0;






/*
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


*/


}
