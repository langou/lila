#include "lila.h"

int lila_dge_qr_larft_connect_w03( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt ){

	double *Aii, *Aij, *T0j, *T0i, *Tmodi;
	int vb, itlo, ithi, jtlo, jthi, jj, ii;

	vb = mt - (i%mt);
	if( vb > n ) vb = n;

	itlo = (i%mt);
	ithi = itlo+vb;
	jtlo = i;
	jthi = jtlo + vb;

	Aii = A + i + i*lda;
	Aij = A + i + (jtlo - itlo)*lda;
	
	T0j = T + (jtlo - itlo)*ldt;
	T0i = T + jtlo*ldt;
	Tmodi = T + itlo + jtlo*ldt;

	for( jj = 0; jj < vb; jj++ ){
		for( ii = 0; ii < itlo; ii++ ){
			T0i[ ii + jj * ldt ] = Aij[  jj + ii * lda  ];
		}
	}

	cblas_dtrmm ( CblasColMajor, CblasRight, CblasLower, CblasNoTrans, CblasUnit, itlo, vb, (+1.0e+00), Aii, lda, T0i, ldt );
	cblas_dgemm ( CblasColMajor, CblasTrans, CblasNoTrans, itlo, vb, m-vb-i, (+1.0e+00), Aij+vb, lda, Aii+vb, lda, (+1.0e+00), T0i, ldt );
	cblas_dtrmm ( CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, itlo, vb, (-1.0e+00), T0j, ldt, T0i, ldt );
	cblas_dtrmm ( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, itlo, vb, (+1.0e+00), Tmodi, ldt, T0i, ldt );

	return 0;

}

/*

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

}

*/
