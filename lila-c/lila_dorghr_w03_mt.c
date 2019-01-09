#include "lila.h"

int lila_dorghr_w03_mt( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork ){

	int vb, jj, j, not_done, info;

	vb = mt - ( i % mt );
	if ( vb > n ) vb = n; 

	not_done = 1;
	jj = 1;
	j = i;

	while( not_done == 1 ){


		//if( jj != 1 ) info = lila_ormhr2_w03_hr ( m, vb, j, mt, A, lda, T, ldt, Q, ldq, work, lwork );
	//lila_dorgh2_3( ml, n, Aii, lda, Tki, ldt, Qii, ldq, work, lwork, S );


	
		if( jj != 1 ) info = lila_ormhr2_w03_hr ( m, vb, j, mt, A, lda, T, ldt, Q, ldq, work, lwork );
		info = lila_dgeqr2_w03_hr               ( m, vb, j, mt, A, lda, T, ldt, Q, ldq, work, lwork );
		info = lila_dlarft_connect_w03          ( m, vb, j, i, mt, A, lda, T, ldt );

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
