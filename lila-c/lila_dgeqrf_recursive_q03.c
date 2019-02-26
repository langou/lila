#include "lila.h"

int lila_dgeqrf_recursive_q03( int panel, int leaf, int nx, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork ){

	int info, nb1, nb2; 

	if ( n < nx ) {
		
		if ( leaf == 0 ){
			info = lila_dgeqrf_q03_mt_l ( panel, m, n, i, mt, A, lda, T, ldt, Q, ldq, work, lwork );
		} 
		if ( leaf == 1 ){
			info = lila_dgeqrf_q03_mt   ( panel, m, n, i, mt, A, lda, T, ldt, Q, ldq, work, lwork );
		} 
		if ( leaf == 2 ){
//			info = lila_dgeqrf_q03_mt_hr( panel, m, n, i, mt, A, lda, T, ldt, Q, ldq, work, lwork );
		}

	} else {

		nb1 = n / 2;
		nb2 = n - nb1;

		info = lila_dgeqrf_recursive_q03( panel, leaf, nx, m, nb1, i, mt, A, lda, T, ldt, Q, ldq, work, lwork );
		info = lila_dgeqrf_recursive_q03( panel, leaf, nx, m, nb2, i+nb1, mt, A, lda, T, ldt, Q, ldq, work, lwork );
		info = lila_dormqrbz_w03( m, nb2, nb1, i, i+nb1, mt, A, lda, Q, ldq, T, ldt, work, lwork );

	}

	return 0;

}