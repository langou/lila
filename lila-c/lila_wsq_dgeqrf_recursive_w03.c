#include "lila.h"

int lila_wsq_dgeqrf_recursive_w03( int panel, int leaf, int nx, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork ){

	int nb1, nb2, lwork1; 

	lwork1 = 0;
	if ( n < nx ) {
		
		if ( leaf == 0 ){
			lwork1 = lila_wsq_dgeqrf_w03_mt_l ( panel, m, n, i, mt, A, lda, T, ldt, Q, ldq, work, lwork ); if( lwork < lwork1 ) lwork = lwork1;
		} 
		if ( leaf == 1 ){
			lwork1 = lila_wsq_dgeqrf_w03_mt   ( panel, m, n, i, mt, A, lda, T, ldt, Q, ldq, work, lwork ); if( lwork < lwork1 ) lwork = lwork1;
		} 
		if ( leaf == 2 ){
			lwork1 = lila_wsq_dgeqrf_w03_mt_hr( panel, m, n, i, mt, A, lda, T, ldt, Q, ldq, work, lwork ); if( lwork < lwork1 ) lwork = lwork1;
		}
		printf(" 1 |  lwork  = %3d,\n",lwork);
		
	} else {

		nb1 = n / 2;
		nb2 = n - nb1;

		lwork1 = lila_wsq_dgeqrf_recursive_w03( panel, leaf, nx, m, nb1, i, mt, A, lda, T, ldt, Q, ldq, work, lwork );     if( lwork < lwork1 ) lwork = lwork1;
		lwork1 = lila_wsq_dormqrf_w03( m, nb2, nb1, i, i+nb1, mt, A, lda, T, ldt, work, lwork );                           if( lwork < lwork1 ) lwork = lwork1;
		lwork1 = lila_wsq_dgeqrf_recursive_w03( panel, leaf, nx, m, nb2, i+nb1, mt, A, lda, T, ldt, Q, ldq, work, lwork ); if( lwork < lwork1 ) lwork = lwork1;
		lwork1 = lila_wsq_dormqrbz_w03( m, nb2, nb1, i, i+nb1, mt, A, lda, Q, ldq, T, ldt, work, lwork );                  if( lwork < lwork1 ) lwork = lwork1;

	}
	printf(" 1 |  lwork  = %3d,\n",lwork);
	return (lwork);

}
