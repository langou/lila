#include "lila.h"

int lila_dgeqrf_levelx_w03( int n_lvl, int i_lvl, int *nb_lvl, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork ){

	int vb, info; 
	int j, k, nb;

	nb = nb_lvl[ i_lvl ];

	if ( nb > n ) vb = n; else vb = nb; 

	k = 0;
	j = i;

	if( i_lvl == n_lvl-1 ){
//	printf("vb (level %d) = %d\n", i_lvl, vb );
	info = lila_dgeqr2_w03b( m, vb, j, mt, A, lda, T, ldt, Q, ldq, work, lwork );
	} else {
//	printf("vb (level %d) = %d\n", i_lvl, vb );
	info = lila_dgeqrf_levelx_w03( n_lvl, i_lvl+1, nb_lvl, m, vb, j, mt, A, lda, T, ldt, Q, ldq, work, lwork );
	}

	k += vb;
	j += vb;
	
	while ( j < i+n ) {

	if ( j+nb > i+n ) vb = i+n-j; else vb = nb; 

	info = lila_dormqrf_w03( m, vb, k, i, j, mt, A, lda, T, ldt, work, lwork );

	if( i_lvl == n_lvl-1 ){
//		printf("vb (level %d) = %d\n", i_lvl, vb );
		info = lila_dgeqr2_w03a( m, vb, j, mt, A, lda, T, ldt, Q, ldq, work, lwork );
	} else {
//		printf("vb (level %d) = %d\n", i_lvl, vb );
		info = lila_dgeqrf_levelx_w03( n_lvl, i_lvl+1, nb_lvl, m, vb, j, mt, A, lda, T, ldt, Q, ldq, work, lwork );
	}

	info = lila_dlarft_connect_w03( m, vb, j, i, mt, A, lda, T, ldt );

//	info = lila_dormqrbz_w00( m, vb, k, i, j, mt, A, lda, Q, ldq, T, ldt, work, lwork );
	info = lila_dormqrbz_w03( m, vb, k, i, j, mt, A, lda, Q, ldq, T, ldt, work, lwork );

 
	k += vb;
	j += vb;

	}

	return 0;

}
