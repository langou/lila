#include "lila.h"

int lila_dgeqrf_levelx( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork ){

	int vb, info, j, k, nb, nb_lvl, n_lvl, i_lvl;
	n_lvl = lila_param[3];
	i_lvl  = lila_param[5];
	nb_lvl = lila_param[6+i_lvl];

	nb = nb_lvl;
	if ( nb > n ) vb = n; else vb = nb; 

	k = 0;
	j = i;

	if( i_lvl == n_lvl-1 ){

		info = lila_dgeqrf_mt( lila_param, m, vb, i, mt, A, lda, T, ldt, Q, ldq, work, lwork );

	} else {

		lila_param[5] += 1; // increasing i_lvl
		info = lila_dgeqrf_levelx( lila_param, m, vb, j, mt, A, lda, T, ldt, Q, ldq, work, lwork );

	}

		k += vb;
		j += vb;
	
	while ( j < i+n ) {

		if ( j+nb > i+n ) vb = i+n-j; else vb = nb; 

		info = lila_dormqrf( lila_param, m, vb, k, i, j, mt, A, lda, T, ldt, work, lwork );

	if( i_lvl == n_lvl-1 ){

		info = lila_dgeqrf_mt( lila_param, m, vb, j, mt, A, lda, T, ldt, Q, ldq, work, lwork );

	} else {

		info = lila_dgeqrf_levelx( lila_param, m, vb, j, mt, A, lda, T, ldt, Q, ldq, work, lwork );

	}

		info = lila_dlarft_connect( lila_param, m, vb, j, i, mt, A, lda, T, ldt );
		info = lila_dormqrbz( lila_param, m, vb, k, i, j, mt, A, lda, Q, ldq, T, ldt, work, lwork );
 
		k += vb;
		j += vb;

	}

	return 0;

}
