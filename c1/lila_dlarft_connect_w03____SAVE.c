#include "lila.h"

int lila_dlarft_connect_w03( int m, int n, int i, int j, int mt, double *A, int lda, double *T, int ldt ){

	double *Aii, *Aij, *T0j, *T0i, *Tmodi;
	int vb, itlo, jtlo, jj, ii;
	int wb;

	vb = mt - (i%mt);
	if( vb > n ) vb = n;

	if ((j/mt) == (i/mt)) wb = i - j; else wb = (mt - (j%mt));

	itlo = ( i % mt );
	jtlo = i;

	if( itlo != 0 ) {

	Aii = A + i + i*lda;
//	Aij = A + i + (jtlo - itlo)*lda;
	Aij = A + i + (i - wb)*lda;
	
	T0j = T + (jtlo - itlo)*ldt;
	T0i = T + jtlo*ldt;

	Tmodi = T + itlo + jtlo*ldt;

	for( jj = 0; jj < vb; jj++ ){
		for( ii = 0; ii < itlo; ii++ ){
			T0i[ ii + jj * ldt ] = Aij[  jj + ii * lda  ];
		}
	}
	
	printf(" m = %2d,   n = %2d,   wb = %2d,   vb = %2d,   i = %2d,   j = %2d,   itlo = %2d,   jtlo = %2d,\n",m,n,wb,vb,i,j,itlo,jtlo);

	cblas_dtrmm ( CblasColMajor, CblasRight, CblasLower, CblasNoTrans, CblasUnit, itlo, vb, (+1.0e+00), Aii, lda, T0i, ldt );
	cblas_dgemm ( CblasColMajor, CblasTrans, CblasNoTrans, itlo, vb, m-vb-i, (+1.0e+00), Aij+vb, lda, Aii+vb, lda, (+1.0e+00), T0i, ldt );
	cblas_dtrmm ( CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, itlo, vb, (-1.0e+00), T0j, ldt, T0i, ldt );
	cblas_dtrmm ( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, itlo, vb, (+1.0e+00), Tmodi, ldt, T0i, ldt );


	}

	return 0;

}
