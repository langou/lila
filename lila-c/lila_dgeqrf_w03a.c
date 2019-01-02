#include "lila.h"

int lila_dgeqrf_w03a( int m, int n, int nb, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork ){

	int vb, j, info, jj;

	j = 0;
	jj = i;
	if ( nb > n ) vb = n; else vb = nb;
	info = lila_dgeqr2_w03a( m, vb, i, mt, A, lda, T, ldt, Q, ldq, work, lwork );
	j += vb;
	jj += vb;
	if ( j+nb > n ) vb = n-j; else vb = nb;
	while( vb!=0 ){
	info = lila_dormqrf_w03( m, vb, j, i, jj, mt, A, lda, T, ldt, work, lwork );
	info = lila_dgeqr2_w03a( m, vb, jj, mt, A, lda, T, ldt, Q, ldq, work, lwork );
	info = lila_dlarft_connect_w03(m, vb, jj, i, mt, A, lda, T, ldt );
	info = lila_dormqrbz_w03( m, vb, j, i, jj, mt, A, lda, Q, ldq, T, ldt, work, lwork );	
	j += vb;
	jj += vb;
	if ( j+nb > n ) vb = n-j; else vb = nb;
	}

	return 0;

}
