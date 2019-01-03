#include "lila.h"

int lila_dgeqrf_w03_mt( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork ){

	double *Tki, *Aii, *Qii;
	int vb, info;
	int j, k, l, ml;

	k = i % mt;

//	Tki = T + k + i*ldt;
	Tki = T + (i%mt) + i*ldt;
	Aii = A + i + i*lda;
	Qii = Q + i + i*ldq;

	ml = m - i;
	vb = mt - k; if ( vb > n ) vb = n;
	
//	info = lila_dgeqr2_w03_l ( m, vb, i, mt, A, lda, T, ldt, Q, ldq, work, lwork );
	info = lila_dgeqr2_w03_3 ( m, vb, i, mt, A, lda, T, ldt, Q, ldq, work, lwork );

	j = i + vb;
	l = vb;
	ml -= vb;

	if( j + mt >= i + n ) vb = n - ( j - i ); else vb = mt;

	while( vb != 0 ){

		info = lila_dormqrf_w03  ( m, vb, l, i, j, mt, A, lda, T, ldt, work, lwork );
		info = lila_dgeqr2_w03_l ( m, vb, j, mt, A, lda, T, ldt, Q, ldq, work, lwork );
		info = lila_dormqrbz_w03 ( m, vb, l, i, j, mt, A, lda, Q, ldq, T, ldt, work, lwork );	

		j += vb;
		l += vb;
		ml -= vb;

		if( j + mt >= i + n ) vb = n - ( j - i ); else vb = mt;

	}

	return 0;
}
