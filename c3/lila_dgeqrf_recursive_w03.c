#include "lila.h"

int lila_dgeqrf_recursive_w03( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork ){

	int info, nb1, nb2; 

//	if ( n < 20 ) {

//		info = lila_dgeqr2_w03a( m, n, i, mt, A, lda, T, ldt, Q, ldq, work, lwork );
//		info = lila_dgeqr2_w03b( m, n, i, mt, A, lda, T, ldt, Q, ldq, work, lwork );

	if ( n < 2 ) {

//		info = lila_dgeqr1_w03a( m, n, i, mt, A, lda, T, ldt, Q, ldq, work, lwork );
//		info = lila_dgeqr1_w03b( m, n, i, mt, A, lda, T, ldt, Q, ldq, work, lwork );

//		info = lila_dgeqr2_w03a( m, n, i, mt, A, lda, T, ldt, Q, ldq, work, lwork );
		info = lila_dgeqr2_w03b( m, n, i, mt, A, lda, T, ldt, Q, ldq, work, lwork );

	} else {

	nb1 = n / 2;
	nb2 = n - nb1;

	info = lila_dgeqrf_recursive_w03( m, nb1, i, mt, A, lda, T, ldt, Q, ldq, work, lwork );

//	info = lila_dormqrf_w00( m, nb2, nb1, i, i+nb1, mt, A, lda, T, ldt, work, lwork );
	info = lila_dormqrf_w03( m, nb2, nb1, i, i+nb1, mt, A, lda, T, ldt, work, lwork );

	info = lila_dgeqrf_recursive_w03( m, nb2, i+nb1, mt, A, lda, T, ldt, Q, ldq, work, lwork );

	info = lila_dlarft_connect_w03( m, nb2, i+nb1, i, mt, A, lda, T, ldt );

//	info = lila_dormqrbz_w00( m, nb2, nb1, i, i+nb1, mt, A, lda, Q, ldq, T, ldt, work, lwork );
	info = lila_dormqrbz_w03( m, nb2, nb1, i, i+nb1, mt, A, lda, Q, ldq, T, ldt, work, lwork );

	}

	return 0;

}
