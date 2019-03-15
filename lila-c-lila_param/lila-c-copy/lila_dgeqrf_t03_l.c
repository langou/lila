#include "lila.h"

int lila_dgeqrf_t03_l( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *tau, double *work, int lwork ){

	int vb, info, j;

	vb = mt - ( i%mt ); if ( vb > n ) vb = n;
	j = i;

	info = lila_dgeqr2_t03_l( m, vb, j, mt, A, lda, T, ldt, tau, work, lwork );		
	
	j   = i + vb;

	if( j + mt >= i + n ) vb = n - ( j - i ); else vb = mt;

	while( vb != 0 ){

		info = lila_dgeqr2_t03_l( m, vb, j, mt, A, lda, T, ldt, tau, work, lwork );		

		j  += vb;

		if( j + mt >= i + n ) vb = n - ( j - i ); else vb = mt;

	}

	return 0;
}
