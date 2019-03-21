#include "lila.h"

int lila_dgeqrf_v03_levelx( int *lila_param, int n_lvl, int i_lvl, int *nb_lvl, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *work, int lwork ){

	int vb, info, j, k, nb;

	nb = nb_lvl[ i_lvl ];
	if ( nb > n ) vb = n; else vb = nb; 

	k = 0;
	j = i;

	if( i_lvl == n_lvl-1 ){

		info = lila_dgeqrf_v03_mt( lila_param, m, vb, i, mt, A, lda, T, ldt, work, lwork );

	} else {
		info = lila_dgeqrf_v03_levelx( lila_param, n_lvl, i_lvl+1, nb_lvl, m, vb, j, mt, A, lda, T, ldt, work, lwork );
	}

		k += vb;
		j += vb;
	
	while ( j < i+n ) {

		if ( j+nb > i+n ) vb = i+n-j; else vb = nb; 

		info = lila_dormqrf_w03( lila_param, m, vb, k, i, j, mt, A, lda, T, ldt, work, lwork );

	if( i_lvl == n_lvl-1 ){

		info = lila_dgeqrf_v03_mt( lila_param, m, vb, j, mt, A, lda, T, ldt, work, lwork );

	} else {
		info = lila_dgeqrf_v03_levelx( lila_param, n_lvl, i_lvl+1, nb_lvl, m, vb, j, mt, A, lda, T, ldt, work, lwork );
	}

		info = lila_dlarft_connect_w03( lila_param, m, vb, j, i, mt, A, lda, T, ldt );
 
		k += vb;
		j += vb;

	}

	return 0;

}
