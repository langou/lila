#include "lila.h"

int lila_query_dormqrf_w03( int m, int n, int k, int i, int j, int mt, double *A, int lda, double *T, int ldt, double *work, int lwork ){

	int ml, vb, lwork1, not_done, jj;

	vb = mt - (i % mt );
	if ( vb > k ) vb = k;
	
	lwork1 = 0;
	ml       = m - i;
	not_done = 1;
	jj       = 1;

	if ( lwork1 < n*mt ) lwork1 = n*mt;

	while( not_done == 1 ){

	if ( lwork1 < (ml-vb)*n ) lwork1 = (ml-vb)*n;

	if ( jj + vb - 1 == k ) {
			not_done = 0; 
	} else {
		ml  -= vb;
		jj  += vb;
		if ( ( jj + mt - 1 ) <= k ) vb = mt; else vb = k - jj + 1;
	}
	}

	return lwork1;

}
