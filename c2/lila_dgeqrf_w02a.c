#include "lila.h"

int lila_dgeqrf_w02a( int m, int n, int i, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork, int nb ){

	int vb, info, mt, k, j;

	mt = -1;

	k = 0;
	j = i;
	if ( nb > n ) vb = n; else vb = nb;
	info =        lila_dgeqr2_w02a( m, vb, i, mt, A, lda, T, ldt, Q, ldq, work, lwork );
	k += vb;
	j += vb;
	if ( k+nb > n ) vb = n-k; else vb = nb;

	while( vb!=0 ){
	info =        lila_dormqrf_z02( m, vb, k, i, j, mt, A, lda, T, ldt, A, lda, work, lwork );
	info =        lila_dgeqr2_w02a( m, vb, j,       mt, A, lda, T, ldt, Q, ldq, work, lwork );
	info = lila_dlarft_connect_w02( m, vb, j, i,    mt, A, lda, T, ldt );
	info =       lila_dormqrbz_w02( m, vb, k, i, j, mt, A, lda, T, ldt, Q, ldq, work, lwork );	
	k += vb;
	j += vb;
	if ( k+nb > n ) vb = n-k; else vb = nb;
	}

	return 0;

}
