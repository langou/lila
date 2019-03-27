#include "lila.h"

int lila_dgeqrf_appendcols( int *lila_param, int m, int n, int k, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork ){


	int info;
 
	if( k != 0 ) info = lila_dormqrf       ( lila_param, m, n, k, 0, k, mt, A, lda, T, ldt, work, lwork );
	lila_dgeqrf_recursive                  ( lila_param, m, n, k, mt, A, lda, T, ldt, Q, ldq, work, lwork );
	if( k != 0 ) info = lila_dlarft_connect( lila_param, m, n, k, 0, mt, A, lda, T, ldt );
	if( k != 0 ) info = lila_dormqrbz      ( lila_param, m, n, k, 0, k, mt, A, lda, Q, ldq, T, ldt, work, lwork );

	return 0;

}


/*

		tau = work;
		Akk = Aii;
		for( kk = 0; kk < n; kk++){ normv = 1 + cblas_ddot(ml-kk-1,Akk+1,1,Akk+1,1); tau[kk] = 2.0e+00 / normv; Akk=Akk+1+lda; }


*/
