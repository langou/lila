#include "lila.h"

int lila_dgeqrf_w02b( int m, int n, int i, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork, int nb ){

	int *S;
	int vb, info, mt, k, j;

	S = (int *) malloc(n * sizeof(int));

	j = i;
	k = 0;
	if ( nb > n ) vb = n; else vb = nb;
	lila_ormhr_w0b( m, vb, i, j, A, lda, T, ldt, Q, ldq, S );
	lila_dorgh2( m, vb, i, j, mt, A, lda, T, ldt, Q, ldq, work, lwork, S );
	j += vb;
	k += vb;
	if ( j+nb > n ) vb = n-j; else vb = nb;
	
	while( vb!=0 ){
		lila_ormhr_w0b( m, vb, i, j, A, lda, T, ldt, Q, ldq, S );
		lila_dorgh2( m, vb, i, j, mt, A, lda, T, ldt, Q, ldq, work, lwork, S );
		info = lila_dlarft_connect_w02(m, vb, j, i, mt, A, lda, T, ldt );
		j += vb;
		k += vb;
		if ( j+nb > n ) vb = n-j; else vb = nb;
	}
	
	free( S );

	return 0;
}



/*

lila_dgeqrf_w02a

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


*/
