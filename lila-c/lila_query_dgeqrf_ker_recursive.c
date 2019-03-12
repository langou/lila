#include "lila.h"

int lila_query_dgeqrf_ker_recursive( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork ){

	int lwork1, nb1, nb2, nx; 

	nx = lila_param[ 3 ];

	if ( n <= nx ) {

		lwork1 = n+i; if (lwork1 > lwork) lwork = lwork1;
		lwork1 = n*n; if (lwork1 > lwork) lwork = lwork1;
		
	} else {
	
		nb1 = n / 2;
		nb2 = n - nb1;

		lwork1 = lila_query_dgeqrf_ker_recursive( lila_param, m, nb1, i, mt, A, lda, T, ldt, Q, ldq, work, lwork );
		lwork1 = nb1*nb2; if (lwork1 > lwork) lwork = lwork1;
		lwork1 = lila_query_dgeqrf_ker_recursive( lila_param, m, nb2, i+nb1, mt, A, lda, T, ldt, Q, ldq, work, lwork );
		lwork1 = nb1*nb2; if (lwork1 > lwork) lwork = lwork1;
		lwork1 = nb1*nb2; if (lwork1 > lwork) lwork = lwork1;

	}
	
	return lwork;

}
