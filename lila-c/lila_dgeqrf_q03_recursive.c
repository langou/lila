#include "lila.h"

int lila_dgeqrf_q03_recursive( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork ){

	int info, nb1, nb2, nx;
	nx = lila_param[3]; 

	if ( n < nx ) {
		
		info = lila_dgeqrf_q03_mt( lila_param, m, n, i, mt, A, lda, T, ldt, Q, ldq, work, lwork );

	} else {

		nb1 = n / 2;
		nb2 = n - nb1;

		info = lila_dgeqrf_q03_recursive( lila_param, m, nb1, i, mt, A, lda, T, ldt, Q, ldq, work, lwork );
		info = lila_dgeqrf_q03_recursive( lila_param, m, nb2, i+nb1, mt, A, lda, T, ldt, Q, ldq, work, lwork );
		info = lila_dormqrbz_w03( m, nb2, nb1, i, i+nb1, mt, A, lda, Q, ldq, T, ldt, work, lwork );

	}

	return 0;

}
