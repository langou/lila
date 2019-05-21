#include "lila.h"

int lila_query_dgeqrf_w03_recursive( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork ){

	int nb1, nb2, lwork1; 

	lwork1 = 0;
	if ( n <= lila_param[3] ) {
		
		lwork1 = lila_query_dgeqrf_w03_mt( lila_param, m, n, i, mt, A, lda, T, ldt, Q, ldq, work, lwork ); if( lwork < lwork1 ) lwork = lwork1;
		
	} else {

		nb1 = n / 2;
		nb2 = n - nb1;

		lwork1 = lila_query_dgeqrf_w03_recursive( lila_param, m, nb1, i, mt, A, lda, T, ldt, Q, ldq, work, lwork );     if( lwork < lwork1 ) lwork = lwork1;
		lwork1 = lila_query_dormqrf_w03( m, nb2, nb1, i, i+nb1, mt, A, lda, T, ldt, work, lwork );                           if( lwork < lwork1 ) lwork = lwork1;
		lwork1 = lila_query_dgeqrf_w03_recursive( lila_param, m, nb2, i+nb1, mt, A, lda, T, ldt, Q, ldq, work, lwork ); if( lwork < lwork1 ) lwork = lwork1;
		lwork1 = lila_query_dormqrbz_w03( m, nb2, nb1, i, i+nb1, mt, A, lda, Q, ldq, T, ldt, work, lwork );                  if( lwork < lwork1 ) lwork = lwork1;

	}

	return (lwork);

}
