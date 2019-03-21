#include "lila.h"

int lila_dgeqrf_t03( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *work, int lwork ){

	double *Aii, *tau;
	int info, ml, t03;

	t03 = lila_param[5];

	Aii = A + i + i*lda;
	ml  = m - i;
	tau = work+n; 

	// necessary, we need v (and tau) to compute T
	// do all of the qr-factorization right away
  	info = LAPACKE_dgeqrf_work( LAPACK_COL_MAJOR, ml, n, Aii, lda, tau, work, lwork ); 
	
	if ( t03 == 0 ){
	//	info = LAPACKE_dlaset( LAPACK_COL_MAJOR, 'A', mt, n, (0e+00), (0e+00), T, ldt );
		info = lila_dlarft_w03( m, n, i, mt, A, lda, T, ldt, tau );
	} else { 
		info = lila_dgeqrf_t03_l( m, n, i, mt, A, lda, T, ldt, tau, work, lwork );		
	}

	return 0;

}
