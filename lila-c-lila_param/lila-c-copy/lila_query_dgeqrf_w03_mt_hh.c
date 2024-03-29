#include "lila.h"

int lila_query_dgeqrf_w03_mt_hh( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork ){

	int vb, j, l, ml, lwork1;

	lwork1 = 0;

	ml = m - i;
	vb = mt - ( i%mt ); if ( vb > n ) vb = n;

	if( lila_param[2] == 0 ){	
		lwork1 = lila_query_dgeqr2_w03_l ( m, vb, i, mt, A, lda, T, ldt, Q, ldq, work, lwork ); if ( lwork < lwork1 ) lwork = lwork1;
	}
	else if( lila_param[2] == 1 ){
		lwork1 = lila_query_dgeqr2_w03_3 ( m, vb, i, mt, A, lda, T, ldt, Q, ldq, work, lwork ); if ( lwork < lwork1 ) lwork = lwork1;
	}
	else if ( lila_param[2] == 2 ){
		lwork1 = lila_query_dgeqr2_w03_hr( m, vb, i, mt, A, lda, T, ldt, Q, ldq, work, lwork ); if ( lwork < lwork1 ) lwork = lwork1;
	}

	j   = i + vb;
	l   = vb;
	ml -= vb;
	if( j + mt >= i + n ) vb = n - ( j - i ); else vb = mt;

	while( vb != 0 ){

		lwork1 = lila_query_dormqrf_w03  ( m, vb, l, i, j, mt, A, lda, T, ldt, work, lwork ); if ( lwork < lwork1 ) lwork = lwork1;

		if( lila_param[2] == 0 ){	
			lwork1 = lila_query_dgeqr2_w03_l ( m, vb, j, mt, A, lda, T, ldt, Q, ldq, work, lwork ); if ( lwork < lwork1 ) lwork = lwork1;
		}
		else if( lila_param[2] == 1 ){
			lwork1 = lila_query_dgeqr2_w03_3 ( m, vb, j, mt, A, lda, T, ldt, Q, ldq, work, lwork ); if ( lwork < lwork1 ) lwork = lwork1;
		}
		else if( lila_param[2] == 2 ){
			lwork1 = lila_query_dgeqr2_w03_hr( m, vb, j, mt, A, lda, T, ldt, Q, ldq, work, lwork ); if ( lwork < lwork1 ) lwork = lwork1;
		}

		lwork1 = lila_query_dormqrbz_w03 ( m, vb, l, i, j, mt, A, lda, Q, ldq, T, ldt, work, lwork ); if ( lwork < lwork1 ) lwork = lwork1;	

		j  += vb;
		l  += vb;
		ml -= vb;

		if( j + mt >= i + n ) vb = n - ( j - i ); else vb = mt;

	}

	return lwork;
}
