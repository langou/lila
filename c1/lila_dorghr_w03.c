#include "lila.h"

int lila_dorghr_w03( int m, int n, int i, int j, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, int *S ){

//  This script has not been set-up for anything. It was made for demo_03 but I never came back to implement it

	if ( nb > n ) vb = n; else vb = nb;
		lila_ormhr_w03( m, vb, i, j, mt, A, lda, T, ldt, Q, ldq, S );
		lila_dorgh2_w03( m, vb, i, j, mt, A, lda, T, ldt, Q, ldq, work, lwork, S );
		j += vb;
	if ( j+nb > n ) vb = n-j; else vb = nb;
	while( vb!=0 ){
		lila_ormhr_w03( m, vb, i, j, mt, A, lda, T, ldt, Q, ldq, S );
		lila_dorgh2_w03( m, vb, i, j, mt, A, lda, T, ldt, Q, ldq, work, lwork, S );
		info = lila_dlarft_connect_w03(m, vb, j, i, mt, A, lda, T, ldt );
		j += vb;
	if ( j+nb > n ) vb = n-j; else vb = nb;
	}


	return 0;

}
