#include "lila.h"

int lila_dgeqrf_w03_mt_l( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork ){

	int vrtq;
	vrtq = lila_param[4];

	// vrtq == 0  ==> VRTQ
	if( vrtq == 0 ){

		double *tau=NULL, *Aii, *Qii;
		int ml, info; 

		tau = work + n;

		Aii = A + i*lda + i;
		Qii = Q + i*ldq + i;

		ml = m - i;

	  	info = LAPACKE_dgeqrf_work( LAPACK_COL_MAJOR, ml, n, Aii, lda, tau, work, lwork );
		info = lila_dlarft_w03( m, n, i, mt, A, lda, T, ldt, tau);
		info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml, n, Aii, lda, Qii, ldq );
		info = LAPACKE_dorgqr_work( LAPACK_COL_MAJOR, ml, n, n, Qii, ldq, tau, work, lwork );
	
	}

	// vrtq == 1  ==> VRT
	if( vrtq == 1 ){

		double *tau=NULL, *Aii;
		int ml, info; 

		tau = work + n;
		Aii = A + i*lda + i;

		ml = m - i;

	  	info = LAPACKE_dgeqrf_work( LAPACK_COL_MAJOR, ml, n, Aii, lda, tau, work, lwork );
		info = lila_dlarft_w03( m, n, i, mt, A, lda, T, ldt, tau);

	}

	// vrtq == 2  ==> Q
	if( vrtq == 2 ){

		double *tau=NULL, *Aii, *Qii;
		int ml, info, lwork1; 

		ml     = m - i;
		tau    = work;
		work   = work + n;
		lwork1 = lwork-n;

		Aii = A + i*lda + i;
		Qii = Q + i*ldq + i;

		info = lila_dT2tau_w03( m, n, i, mt, T, ldt, tau);
		info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml, n, Aii, lda, Qii, ldq );
		info = LAPACKE_dorgqr_work( LAPACK_COL_MAJOR, ml, n, n, Qii, ldq, tau, work, lwork1 );

	}

	return 0;

}
