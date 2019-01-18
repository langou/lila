#include "lila.h"

int lila_dlarft_connect_w03( int m, int n, int i, int j, int mt, double *A, int lda, double *T, int ldt ){

	double *Aii, *Aij, *Tjj, *Tji, *Tmodi;
	int vb, itlo, jtlo, jj, ii, wb;

	itlo = ( i % mt );
	vb = mt - (i%mt);

	if( vb > n ) vb = n;
	if( (j/mt) == (i/mt)) wb = i - j; else wb = (i%mt);
	if( wb == mt ) wb = 1; 

	Aii = A + i + i*lda;
	Tmodi = T + itlo + i*ldt;

	Aij = A + i + (i-wb)*lda;
	Tji = T + (itlo-wb) + i*ldt;
	Tjj = T + (itlo-wb) + (i-wb)*ldt;

	for( jj = 0; jj < vb; jj++ ){
		for( ii = 0; ii < wb; ii++ ){
			Tji[ ii + jj * ldt ] = Aij[  jj + ii * lda  ];
		}
	}
	
	cblas_dtrmm ( CblasColMajor, CblasRight, CblasLower, CblasNoTrans, CblasUnit, wb, vb, (+1.0e+00), Aii, lda, Tji, ldt );
	cblas_dgemm ( CblasColMajor, CblasTrans, CblasNoTrans, wb, vb, m-vb-i, (+1.0e+00), Aij+vb, lda, Aii+vb, lda, (+1.0e+00), Tji, ldt );
	cblas_dtrmm ( CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, wb, vb, (-1.0e+00), Tjj, ldt, Tji, ldt );
	cblas_dtrmm ( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, wb, vb, (+1.0e+00), Tmodi, ldt, Tji, ldt );
	
	return 0;

}
