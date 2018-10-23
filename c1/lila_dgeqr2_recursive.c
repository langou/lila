#include "lila.h"

int lila_dgeqr2_recursive( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *TTT, int llldddttt, double *Q, int ldq, double *work, int lwork ){

	int info, nb1, nb2; 


	if ( n < 2 ) {

		printf(" i = %d, (leaf) n = %d,\n", i,n);

//		info = lila_dge_qr_wq_w03( m, n, i, mt, A, lda, T, ldt, TTT, llldddttt, Q, ldq, work, lwork );
		info = lila_dgeqr2_w02b( m, n, i, mt, A, lda, T, ldt, Q, ldq, work, lwork ); // This function doesn't have TTT llldddttt in the interface

	} else {

	nb1 = n / 2;
	nb2 = n - nb1;
	printf(" i = %d, nb1 = %d, nb2 = %d, n = %d,\n", i,nb1,nb2,n);

//	info = lila_dge_qr_wq_w03( m, nb1, i, mt, A, lda, T, ldt, TTT, llldddttt, Q, ldq, work, lwork );
//	info = lila_dgeqr2_w02b( m, nb1, i, mt, A, lda, T, ldt, Q, ldq, work, lwork ); // This function doesn't have TTT llldddttt in the interface
	info = lila_dgeqr2_recursive( m, nb1, i, mt, A, lda, T, ldt, TTT, llldddttt, Q, ldq, work, lwork );

	info = lila_dge_qr_ormqrf_w03( m, nb2, nb1, i, i+nb1, mt, A, lda, T, ldt, TTT, llldddttt, work, lwork );

	info = lila_dgeqr2_recursive( m, nb2, i+nb1, mt, A, lda, T, ldt, TTT, llldddttt, Q, ldq, work, lwork );

	lila_dge_qr_larft_connect_w02( m, nb2, i+nb1, mt, A, lda, T, ldt );
	lila_dge_qr_larft_connect_w03( m, nb2, i+nb1, mt, A, lda, TTT, llldddttt );

	info = lila_dge_qr_ormqrbz_w00( m, nb2, nb1, i, i+nb1, mt, A, lda, Q, ldq, T, ldt, work, lwork );
//	This one will need to wait until I get it working completely.
//	info = lila_dge_qr_ormqrbz_w03( m, nb2, nb1, i, i+nb1, mt, A, lda, Q, ldq, TTT, llldddttt, work, lwork );
	}

	return 0;

}
