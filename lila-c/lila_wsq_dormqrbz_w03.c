#include "lila.h"

int lila_wsq_dormqrbz_w03( int m, int n, int k, int i, int j, int mt, double *A, int lda, double *Q, int ldq, double *T, int ldt, double *work, int lwork ){

	int vb, jj, lwork1, not_done, ml;

	vb = ((i+k) % mt); if (vb == 0) vb = mt; if ( vb > k ) vb = k;
	lwork1 = 0;
	
	if( lwork1 < n*k ) lwork1 = n*k;

	ml = m - (i+k-vb);
	not_done = 1;
	jj = 1;

	while ( not_done == 1 ){
		if ( jj + vb - 1 == k ){ 
			not_done = 0; 
		} else{

			if( lwork1 < n*vb ) lwork1 = n*vb;
			if( lwork1 < (ml-vb)*k ) lwork1 = (ml-vb)*k;
			//if( lwork1 < ml-vb ) lwork1 = ml-vb;

			jj += vb;
			if ( ( jj + mt - 1 ) <= k ) vb = mt; else { vb = k - jj + 1;  }
			ml += vb;
		}
	}

	return lwork1;

}
