#include "lila.h"

int lila_query_dgeqrf_w03_mt( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork ){

	int lwork1;

		if ( lila_param[1] == 0 ){
			lwork1 = lila_query_dgeqrf_w03_mt_l  ( m, n, i, mt, A, lda, T, ldt, Q, ldq, work, lwork ); if( lwork < lwork1 ) lwork = lwork1;
		} 
		if ( lila_param[1] == 1 ){
			lwork1 = lila_query_dgeqrf_w03_mt_hh ( lila_param, m, n, i, mt, A, lda, T, ldt, Q, ldq, work, lwork ); if( lwork < lwork1 ) lwork = lwork1;
		} 
		if ( lila_param[1] == 2 ){
			lwork1 = lila_query_dgeqrf_w03_mt_hr ( m, n, i, mt, A, lda, T, ldt, Q, ldq, work, lwork ); if( lwork < lwork1 ) lwork = lwork1;
		}

	return lwork;
}
