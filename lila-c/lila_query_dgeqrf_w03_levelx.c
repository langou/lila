#include "lila.h"

int lila_query_dgeqrf_w03_levelx( int *lila_param, int n_lvl, int i_lvl, int *nb_lvl, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork ){

	int vb, j, k, nb, lwork1;

	lwork  = 0;
	lwork1 = 0;

	nb = nb_lvl[ i_lvl ];
	if ( nb > n ) vb = n; else vb = nb; 

	k = 0;
	j = i;

	if( i_lvl == n_lvl-1 ){

		lwork1 = lila_query_dgeqrf_w03_mt( lila_param, m, n, i, mt, A, lda, T, ldt, Q, ldq, work, lwork ); if ( lwork1 > lwork ) lwork = lwork1;

	} else {

		lwork1 = lila_query_dgeqrf_w03_levelx( lila_param, n_lvl, i_lvl+1, nb_lvl, m, vb, j, mt, A, lda, T, ldt, Q, ldq, work, lwork ); if ( lwork1 > lwork ) lwork = lwork1;
	}

		//printf(" 1 |  lwork  = %3d,\n",lwork);

		k += vb;
		j += vb;
	
	while ( j < i+n ) {

		if ( j+nb > i+n ) vb = i+n-j; else vb = nb; 

		lwork1 = lila_query_dormqrf_w03( m, vb, k, i, j, mt, A, lda, T, ldt, work, lwork ); if ( lwork1 > lwork ) lwork = lwork1;

	if( i_lvl == n_lvl-1 ){

		lwork1 = lila_query_dgeqrf_w03_mt( lila_param, m, n, i, mt, A, lda, T, ldt, Q, ldq, work, lwork ); if ( lwork1 > lwork ) lwork = lwork1;

	} else {

		lwork1 = lila_query_dgeqrf_w03_levelx( lila_param, n_lvl, i_lvl+1, nb_lvl, m, vb, j, mt, A, lda, T, ldt, Q, ldq, work, lwork ); if ( lwork1 > lwork ) lwork = lwork1;
	}

		lwork1 = lila_query_dormqrbz_w03( m, vb, k, i, j, mt, A, lda, Q, ldq, T, ldt, work, lwork ); if ( lwork1 > lwork ) lwork = lwork1;
 
		k += vb;
		j += vb;

	}

	//printf(" 1-|  lwork  = %3d,\n",lwork);
	return lwork;

}


