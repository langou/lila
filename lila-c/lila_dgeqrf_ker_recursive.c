#include "lila.h"

int lila_dgeqrf_ker_recursive( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork ){

	int info, nb1, nb2, ml, nx; 
	double *tau=NULL, *Aii, *Qii;

	nx = lila_param[ 3 ];
	tau = T+i;

	Aii = A + i + i*lda;
	Qii = Q + i + i*ldq;

	ml = m - i;

	if ( n <= nx ) {

  		info = LAPACKE_dgeqrf_work( LAPACK_COL_MAJOR, ml, n, Aii, lda, tau, work, lwork );
		info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml, n, Aii, lda, Qii, ldq );
		info = LAPACKE_dorgqr_work( LAPACK_COL_MAJOR, ml, n, n, Qii, ldq, tau, work, lwork );
		
	} else {
	
		nb1 = n / 2;
		nb2 = n - nb1;

		info = lila_dgeqrf_ker_recursive( lila_param, m, nb1, i, mt, A, lda, T, ldt, Q, ldq, work, lwork );
		info = LAPACKE_dormqr_work( LAPACK_COL_MAJOR, 'L', 'T', ml, nb2, nb1, Aii, lda, tau, Aii+nb1*lda, ldq, work, lwork );
		info = lila_dgeqrf_ker_recursive( lila_param, m, nb2, i+nb1, mt, A, lda, T, ldt, Q, ldq, work, lwork );
		
		info = LAPACKE_dlaset( LAPACK_COL_MAJOR, 'A', nb1, nb2, (0e+00), (0e+00), Qii+nb1*ldq, ldq );
		info = LAPACKE_dormqr_work( LAPACK_COL_MAJOR, 'L', 'N', ml, nb2, nb1, Aii, lda, tau, Qii+nb1*ldq, ldq, work, lwork );

	}
	
	return 0;

}
