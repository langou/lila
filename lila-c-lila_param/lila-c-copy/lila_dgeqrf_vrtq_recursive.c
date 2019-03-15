#include "lila.h"

int lila_dgeqrf_vrtq_recursive( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork ){

	int vrtq, info;
	vrtq = lila_param[4];

	if         ( vrtq == 0 ){
		info = lila_dgeqrf_w03_recursive( lila_param, m, n, i, mt, A, lda, T, ldt, Q, ldq, work, lwork );
	} else if ( vrtq == 1 ) {
		info = lila_dgeqrf_v03_recursive( lila_param, m, n, i, mt, A, lda, T, ldt, work, lwork );
	} else if ( vrtq == 2 ) {
		info = lila_dgeqrf_v03_recursive( lila_param, m, n, i, mt, A, lda, T, ldt, work, lwork );
		info = lila_dgeqrf_q03_recursive( lila_param, m, n, i, mt, A, lda, T, ldt, Q, ldq, work, lwork );
	} else if ( vrtq == 3 ) {
		info = lila_dgeqrf_t03          ( lila_param, m, n, i, mt, A, lda, T, ldt, work, lwork );
	}

	return 0;

}
