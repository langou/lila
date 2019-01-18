#include "lila.h"

int lila_dlarft_w03_b( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt ){

	double *Tki, *Aii;
	int vb, k, j, info;

	k = i % mt;

	Aii = A + i + i*lda;
	Tki = T + k + i*ldt;

	vb = mt - k; if ( vb > n ) vb = n;

//	this could be rewritten to use less flops
//	info = LAPACKE_dlaset( LAPACK_COL_MAJOR, 'L', vb-1, vb-1, (0e+00), (0e+00), Tki+1, ldt );
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'U', vb, vb, Aii, lda, Tki, ldt ); 
	cblas_dtrsm( CblasColMajor, CblasRight, CblasLower, CblasTrans, CblasUnit, vb, vb, 1.0e+00, Aii, lda, Tki, ldt );

	j = i + vb;
	
	Aii += vb*(1+lda);
	Tki = T + j*ldt;

	if( j+mt >= i+n ) vb = n-(j-i); else vb = mt;

	while( vb != 0 ){

//		info = LAPACKE_dlaset( LAPACK_COL_MAJOR, 'L', vb-1, vb-1, (0e+00), (0e+00), Tki+1, ldt );
		info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'U', vb, vb, Aii, lda, Tki, ldt ); 
		cblas_dtrsm( CblasColMajor, CblasRight, CblasLower, CblasTrans, CblasUnit, vb, vb, 1.0e+00, Aii, lda, Tki, ldt );

		j += vb;

		Tki += ( vb*ldt );
		Aii += vb*(1+lda);

		if( j+mt >= i+n ) vb = n-(j-i); else vb = mt;
	}

	return 0;
}   
