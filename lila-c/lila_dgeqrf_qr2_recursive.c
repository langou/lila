#include "lila.h"

int lila_dgeqrf_qr2_recursive( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork ){


	int info, nb1, nb2; 

	if ( n <= lila_param[ 3 ] ) {
		
		info = lila_dgeqrf_qr2 ( lila_param, m, n, i, mt, A, lda, T, ldt, Q, ldq, work, lwork );


	} else {

		nb1 = n / 2;
		nb2 = n - nb1;

		info = lila_dgeqrf_qr2_recursive2( lila_param, m, n, i, mt, nb1, nb2, A, lda, T, ldt, Q, ldq, work, lwork );
	}

	return 0;

}
