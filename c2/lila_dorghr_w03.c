#include "lila.h"

int lila_dorghr_w03( int m, int n, int i, int j, int l, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork, int *S ){

	int vb, jj, not_done, info;

	vb = mt - ( j % mt );
	if ( vb > n ) vb = n; 

	not_done = 1;
	jj = 1;

	while( not_done == 1 ){


		lila_ormhr_w03( m, vb, i, j, l, mt, A, lda, T, ldt, Q, ldq, work, lwork, S );
		lila_dorgh2_w03( m, vb, i, j, l, mt, A, lda, T, ldt, Q, ldq, work, lwork, S );
		info = lila_dlarft_connect_w03(m, vb, j, i, mt, A, lda, T, ldt );

		if ( jj + vb - 1 == n ) {
		
			not_done = 0;
		
		} else if ( jj + vb - 1 > n ) {

			printf("defensive programming, we should never have been there, abort\n"); return 0;
	
		} else {

			jj += vb;
			j += vb;

			if( jj + mt - 1 <= n ) vb = mt; else vb = n - jj + 1;
		
		}


	}


	return 0;

}
