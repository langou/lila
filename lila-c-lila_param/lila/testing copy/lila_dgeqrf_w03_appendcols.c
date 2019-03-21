#include "lila.h"

int lila_dgeqrf_w03_appendcols( int *lila_param, int m, int k, int n, int mt, double *A, int lda, double *Q, int ldq, double *T, int ldt, double *work, int lwork ){

	int info; 

	if( k != 0 ) info = lila_dormqrf_w03       ( m, n, k, 0, k, mt, A, lda, T, ldt, work, lwork );
	lila_dgeqrf_w03_recursive                  ( lila_param, m, n, k, mt, A, lda, T, ldt, Q, ldq, work, lwork );
	if( k != 0 ) info = lila_dlarft_connect_w03( m, n, k, 0, mt, A, lda, T, ldt );
	if( k != 0 ) info = lila_dormqrbz_w03      ( m, n, k, 0, k, mt, A, lda, Q, ldq, T, ldt, work, lwork );
	
	return 0;

}
