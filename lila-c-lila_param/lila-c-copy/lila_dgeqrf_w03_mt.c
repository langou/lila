#include "lila.h"

int lila_dgeqrf_w03_mt( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork ){

	int info, leaf;
	leaf = lila_param[1];

	if ( leaf == 0 ){
		info = lila_dgeqrf_w03_mt_l  ( m, n, i, mt, A, lda, T, ldt, Q, ldq, work, lwork );
	} 
	if ( leaf == 1 ){
		info = lila_dgeqrf_w03_mt_hh ( lila_param, m, n, i, mt, A, lda, T, ldt, Q, ldq, work, lwork );
	} 
	if ( leaf == 2 ){
		info = lila_dgeqrf_w03_mt_hr ( m, n, i, mt, A, lda, T, ldt, Q, ldq, work, lwork );
	}

	return 0;
}
