#include "lila.h"

int lila_dlarft_w03( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *tau ){

	double *Tij, *Aii;
	int vb, info;
	int jtlo, ml;

	jtlo = i;

	Tij = T + (i%mt) + i*ldt;
	Aii = A + i + i*lda;

	ml = m - i;

	vb = mt - (i%mt);
	if ( vb > n ) vb = n;

	info = LAPACKE_dlarft_work( LAPACK_COL_MAJOR, 'F', 'C', ml, vb, Aii, lda, tau, Tij, ldt);

	jtlo += vb;

	Aii += ( vb*(1+lda) ) ;
	Tij = T + jtlo*ldt;
	tau += vb;

	ml -= vb;

	if( jtlo+mt >= i+n ) vb = n - ( jtlo - i ); else vb = mt;

	while( vb != 0 ){

		info = LAPACKE_dlarft_work( LAPACK_COL_MAJOR, 'F', 'C', ml, vb, Aii, lda, tau, Tij, ldt);
		
		jtlo += vb;

		Aii += ( vb*(1+lda) ) ;
		Tij += ( vb*ldt );
		tau += vb;

		ml -= vb;

		if( jtlo+mt >= i+n ) vb = n - ( jtlo - i ); else vb = mt;

	}

	return 0;
}   
