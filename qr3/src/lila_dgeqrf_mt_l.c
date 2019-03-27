#include "lila.h"

int lila_dgeqrf_mt_l( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork ){

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
		info = lila_dlarft        ( lila_param, m, n, i, mt, A, lda, T, ldt, tau);
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
		info = lila_dlarft        ( lila_param, m, n, i, mt, A, lda, T, ldt, tau);

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

		info = lila_dT2tau        ( lila_param, m, n, i, mt, T, ldt, tau);
		info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml, n, Aii, lda, Qii, ldq );
		info = LAPACKE_dorgqr_work( LAPACK_COL_MAJOR, ml, n, n, Qii, ldq, tau, work, lwork1 );

	}

	// vrtq == 3  ==> T
	if( vrtq == 3 ){

		double *tau=NULL, *Aii, *Akk, normv2;
		int ml, info, k; 

		tau = work + n;
		Aii = A + i*lda + i;

		ml  = m - i;
		Akk = Aii;
		for( k = 0; k < n; k++){ normv2=1+cblas_ddot(ml-k-1,Akk+1,1,Akk+1,1); tau[k] = 2.0e+00 / normv2; Akk=Akk+1+lda; }

		ml = m - i;

		info = lila_dlarft( lila_param, m, n, i, mt, A, lda, T, ldt, tau);

	}

	return 0;

}
