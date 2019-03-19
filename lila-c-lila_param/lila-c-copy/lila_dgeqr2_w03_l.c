#include "lila.h"

int lila_dgeqr2_w03_l( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork ){

	int vrtq;
	vrtq = lila_param[4];
	
	// vrtq == 0  ==> VRTQ
	if( vrtq == 0 ){

		double *tau=NULL, *Aii, *Qii, *Tki;
		int k, ml, info, lwork1; 

		tau = work;
		work = work + n;
		lwork1 = lwork-n;

		ml = m-i;
		k  = i%mt;

		Aii = A + i + i*lda;
		Qii = Q + i + i*ldq;
		Tki = T + k + i*ldt;

	  	info = LAPACKE_dgeqrf_work( LAPACK_COL_MAJOR, ml, n, Aii, lda, tau, work, lwork1 ); 
		info = LAPACKE_dlarft_work( LAPACK_COL_MAJOR, 'F', 'C', ml, n, Aii, lda, tau, Tki, ldt);
		info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml, n, Aii, lda, Qii, ldq ); 
		info = LAPACKE_dorgqr_work( LAPACK_COL_MAJOR, ml, n, n, Qii, ldq, tau, work, lwork1 );

	}

	// vrtq == 1  ==> VRT
	if( vrtq == 1 ){

		double *tau=NULL, *Aii, *Tki;
		int k, ml, info, lwork1; 

		tau = work;
		work = work + n;
		lwork1 = lwork-n;

		ml = m-i;
		k  = i%mt;

		Aii = A + i + i*lda;
		Tki = T + k + i*ldt;

	  	info = LAPACKE_dgeqrf_work( LAPACK_COL_MAJOR, ml, n, Aii, lda, tau, work, lwork1 ); 
		info = LAPACKE_dlarft_work( LAPACK_COL_MAJOR, 'F', 'C', ml, n, Aii, lda, tau, Tki, ldt);

	}

	// vrtq == 2  ==> Q
	if( vrtq == 2 ){
	
		double *tau=NULL, *Aii, *Qii, *Tki;
		int k, ml, info, lwork1; 

		tau = work;
		work = work + n;
		lwork1 = lwork-n;

		ml = m-i;
		k  = i%mt;

		Aii = A + i + i*lda;
		Qii = Q + i + i*ldq;
		Tki = T + k + i*ldt;

		for (k=0;k<n;k++) tau[k]=Tki[k+k*ldt];
		info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml, n, Aii, lda, Qii, ldq ); 
		info = LAPACKE_dorgqr_work( LAPACK_COL_MAJOR, ml, n, n, Qii, ldq, tau, work, lwork1 );

	}

	return 0;

}
