#include "lila.h"

int lila_dlarft_connect_w02( int m, int n, int i, int j, int mt, double *A, int lda, double *T, int ldt ){

	double *Aii, *Aij, *Tii, *Tji, *Tjj;
	int ii, jj;

	Aii = A + i + i*lda;
	Tii = T + i + i*ldt;

	Aij = A + i + j*lda;
	Tji = T + j + i*ldt;
	Tjj = T + j + j*ldt;

	for( jj = 0; jj < n; jj++ ){
	for( ii = 0; ii < i-j; ii++ ){
		Tji[ ii + jj * ldt ] = Aij[  jj + ii * lda  ];
	}}

	if( i > j ) cblas_dtrmm ( CblasColMajor, CblasRight, CblasLower, CblasNoTrans, CblasUnit, i-j, n, (+1.0e+00), Aii, lda, Tji, ldt );
	if( i > j ) cblas_dgemm ( CblasColMajor, CblasTrans, CblasNoTrans, i-j, n, m-n-i, (+1.0e+00), Aij+n, lda, Aii+n, lda, (+1.0e+00), Tji, ldt );
	if( i > j ) cblas_dtrmm ( CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, i-j, n, (-1.0e+00), Tjj, ldt, Tji, ldt );
	if( i > j ) cblas_dtrmm ( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, i-j, n, (+1.0e+00), Tii, ldt, Tji, ldt );

	return 0;

}
