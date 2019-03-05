#include "lila.h"

int lila_wsq_dgeqrf_w03_mt_hr( int panel, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork ){

	double *S;
	int vb, j, lwork1;
	
	S = work;	

	lwork1 = 0;
	j  = i;
	vb = mt - ( i%mt ); if ( vb > n ) vb = n;

	while( vb != 0 ){

		lwork1 = lila_wsq_ormhr2_w03_hr( m, vb, i, j, -1, mt, A, lda, T, ldt, Q, ldq, work, lwork, S ); if( lwork < lwork1 ) lwork = lwork1;
		j += vb;
		if( j + mt >= i + n ) vb = n - ( j - i ); else vb = mt;

	}


	return lwork;

}
