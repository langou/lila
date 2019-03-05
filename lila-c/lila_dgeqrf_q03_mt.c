#include "lila.h"

int lila_dgeqrf_q03_mt( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork ){

	int info;

	if ( lila_param[1] == 0 ){
		info = lila_dgeqrf_q03_mt_l ( m, n, i, mt, A, lda, T, ldt, Q, ldq, work, lwork );
	} 
	else if ( lila_param[1] == 1 ){
		info = lila_dgeqrf_q03_mt_hh( lila_param, m, n, i, mt, A, lda, T, ldt, Q, ldq, work, lwork );
	} 
	else if ( lila_param[1] == 2 ){
//		info = lila_dgeqrf_q03_mt_hr( panel, m, vb, i, mt, A, lda, T, ldt, Q, ldq, work, lwork );
	}

	return 0;
}
