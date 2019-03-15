#include "lila.h"

int lila_dgeqrf_vrtq_levelx( int *lila_param, int n_lvl, int i_lvl, int *nb_lvl, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork ){

	int vrtq, info;
	vrtq = lila_param[4];

	if         ( vrtq == 0 ){
		info = lila_dgeqrf_w03_levelx( lila_param, n_lvl, i_lvl, nb_lvl, m, n, i, mt, A, lda, T, ldt, Q, ldq, work, lwork );
	} else if ( vrtq == 1 ) {
		info = lila_dgeqrf_v03_levelx( lila_param, n_lvl, i_lvl, nb_lvl, m, n, i, mt, A, lda, T, ldt, work, lwork );
	} else if ( vrtq == 2 ) {
		info = lila_dgeqrf_v03_levelx( lila_param, n_lvl, i_lvl, nb_lvl, m, n, i, mt, A, lda, T, ldt, work, lwork );
		info = lila_dgeqrf_q03_levelx( lila_param, n_lvl, i_lvl, nb_lvl, m, n, i, mt, A, lda, T, ldt, Q, ldq, work, lwork );
	}

	return 0;


}
