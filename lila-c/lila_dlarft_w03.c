#include "lila.h"

int lila_dlarft_w03( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *tau ){

	double *Tki, *Aii, *T0j, *Ajj;
	int vb, info, k, j, ml;

	k = i % mt;

	Tki = T + k + i*ldt;
	Aii = A + i + i*lda;

	ml = m - i;
	vb = mt - k; if ( vb > n ) vb = n;

	info = LAPACKE_dlarft_work( LAPACK_COL_MAJOR, 'F', 'C', ml, vb, Aii, lda, tau, Tki, ldt);

	j = i + vb;
	ml -= vb;

	Ajj = A + j + j*lda;
	T0j = T     + j*ldt;
	tau += vb;

	if( j + mt >= i + n ) vb = n - ( j - i ); else vb = mt;

	while( vb != 0 ){

		info = LAPACKE_dlarft_work( LAPACK_COL_MAJOR, 'F', 'C', ml, vb, Ajj, lda, tau, T0j, ldt);
		
		j += vb;
		ml -= vb;

		Ajj += ( vb + vb * lda ) ;
		T0j += ( vb*ldt );
		tau += vb;

		if( j + mt >= i + n ) vb = n - ( j - i ); else vb = mt;

	}

	return 0;
}   
