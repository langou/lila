#include "lila.h"

int lila_dgeqrf_q03_mt_hh( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork ){

	int vb, info, j, l, ml;

	ml = m - i;
	vb = mt - ( i%mt ); if ( vb > n ) vb = n;

	if( lila_param[2] == 0 ){	
		info = lila_dgeqr2_q03_l ( m, vb, i, mt, A, lda, T, ldt, Q, ldq, work, lwork );
	}
	else if( lila_param[2] == 1 ){
		info = lila_dgeqr2_q03_3 ( m, vb, i, mt, A, lda, T, ldt, Q, ldq, work, lwork );
	}
	else if( lila_param[2] == 2){
//		info = lila_dgeqr2_w03_hr( m, vb, i, mt, A, lda, T, ldt, Q, ldq, work, lwork );
	}

	j   = i + vb;
	l   = vb;
	ml -= vb;

	if( j + mt >= i + n ) vb = n - ( j - i ); else vb = mt;

	while( vb != 0 ){

		if( lila_param[2] == 0 ){	
			info = lila_dgeqr2_q03_l ( m, vb, j, mt, A, lda, T, ldt, Q, ldq, work, lwork );
		}
		else if( lila_param[2] == 1 ){
			info = lila_dgeqr2_q03_3 ( m, vb, j, mt, A, lda, T, ldt, Q, ldq, work, lwork );
		}
		else if( lila_param[2] == 2){
//			info = lila_dgeqr2_w03_hr( m, vb, j, mt, A, lda, T, ldt, Q, ldq, work, lwork );
		}

		info = lila_dormqrbz_w03 ( m, vb, l, i, j, mt, A, lda, Q, ldq, T, ldt, work, lwork );	

		j  += vb;
		l  += vb;
		ml -= vb;

		if( j + mt >= i + n ) vb = n - ( j - i ); else vb = mt;

	}

	return 0;
}
