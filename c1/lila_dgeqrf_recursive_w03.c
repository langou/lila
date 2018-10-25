#include "lila.h"

int lila_dgeqrf_recursive_w03( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork ){

	int info, nb1, nb2; 

	if ( n < 2 ) {

		printf(" i = %d, (leaf) n = %d,\n", i,n);

//		info = lila_dgeqr1_w03a( m, n, i, mt, A, lda, T, ldt, Q, ldq, work, lwork );
//		info = lila_dgeqr1_w03b( m, n, i, mt, A, lda, T, ldt, Q, ldq, work, lwork );
		info = lila_dgeqr2_w03b( m, n, i, mt, A, lda, T, ldt, Q, ldq, work, lwork );

	} else {

	nb1 = n / 2;
	nb2 = n - nb1;
	printf(" i = %d, nb1 = %d, nb2 = %d, n = %d,\n", i,nb1,nb2,n);

	info = lila_dgeqrf_recursive_w03( m, nb1, i, mt, A, lda, T, ldt, Q, ldq, work, lwork );

//	info = lila_dge_qr_ormqrf_w00( m, nb2, nb1, i, i+nb1, mt, A, lda, T, ldt, work, lwork );
	info = lila_dge_qr_ormqrf_w03( m, nb2, nb1, i, i+nb1, mt, A, lda, T, ldt, work, lwork );

	info = lila_dgeqrf_recursive_w03( m, nb2, i+nb1, mt, A, lda, T, ldt, Q, ldq, work, lwork );

	lila_dge_qr_larft_connect_w03( m, nb2, i+nb1, mt, A, lda, T, ldt );

//	info = lila_dge_qr_ormqrbz_w00( m, nb2, nb1, i, i+nb1, mt, A, lda, Q, ldq, T, ldt, work, lwork );
	info = lila_dge_qr_ormqrbz_w03( m, nb2, nb1, i, i+nb1, mt, A, lda, Q, ldq, T, ldt, work, lwork );

	}

	return 0;

}
