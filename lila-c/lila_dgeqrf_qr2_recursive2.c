#include "lila.h"

int lila_dgeqrf_qr2_recursive2( int *lila_param, int m, int n, int i, int mt, int nb1, int nb2, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork ){


	int info; 

	if( lila_param[4] == 0 ){

		info = lila_dgeqrf_qr2_recursive( lila_param, m, nb1, i, mt, A, lda, T, ldt, Q, ldq, work, lwork );
		info = lila_dormqrf_w03( m, nb2, nb1, i, i+nb1, mt, A, lda, T, ldt, work, lwork );
		info = lila_dgeqrf_qr2_recursive( lila_param, m, nb2, i+nb1, mt, A, lda, T, ldt, Q, ldq, work, lwork );
		info = lila_dlarft_connect_w03( m, nb2, i+nb1, i, mt, A, lda, T, ldt );
		info = lila_dormqrbz_w03( m, nb2, nb1, i, i+nb1, mt, A, lda, Q, ldq, T, ldt, work, lwork );
	
	} else if( lila_param[4] == 1){

		info = lila_dgeqrf_qr2_recursive( lila_param, m, nb1, i, mt, A, lda, T, ldt, Q, ldq, work, lwork );
		info = lila_dormqrf_w03( m, nb2, nb1, i, i+nb1, mt, A, lda, T, ldt, work, lwork );
		info = lila_dgeqrf_qr2_recursive( lila_param, m, nb2, i+nb1, mt, A, lda, T, ldt, Q, ldq, work, lwork );
		info = lila_dlarft_connect_w03( m, nb2, i+nb1, i, mt, A, lda, T, ldt );

	} else if( lila_param[4] == 2){

		info = lila_dgeqrf_qr2_recursive( lila_param, m, nb1, i, mt, A, lda, T, ldt, Q, ldq, work, lwork );
		info = lila_dgeqrf_qr2_recursive( lila_param, m, nb2, i+nb1, mt, A, lda, T, ldt, Q, ldq, work, lwork );
		info = lila_dormqrbz_w03( m, nb2, nb1, i, i+nb1, mt, A, lda, Q, ldq, T, ldt, work, lwork );

	}

	return 0;

}
