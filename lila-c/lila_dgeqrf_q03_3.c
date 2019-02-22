#include "lila.h"

int lila_dgeqrf_q03_3( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork ){

	double *Aii, *Qii, *Tki, *tau;
	int vb, info, j, l, ml, k, lwork1;

	tau  = work;
	work = work + n;
	lwork1 = lwork - n;
 
	Aii = A + i      + i*lda;
	Qii = Q + i      + i*ldq;
	Tki = T + (i%mt) + i*ldt;

	ml = m  - i;
	vb = mt - ( i%mt ); if ( vb > n ) vb = n;
	k  = 0;

	for (i=0;i<vb;i++) tau[i]=Tki[(i%mt)+i*ldt];
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml, n, Aii, lda, Qii, ldq ); 
	info = LAPACKE_dorgqr_work( LAPACK_COL_MAJOR, ml, n, n, Qii, ldq, tau, work, lwork1 );

//	info = lila_dormqrbz_w03( m, vb, k, i, i, mt, A, lda, Q, ldq, T, ldt, work, lwork1 );	

	k   += vb;
	j   = i + vb;
	l   = vb;
	ml -= vb;

	if( j + mt >= i + n ) vb = n - ( j - i ); else vb = mt;

	while( vb != 0 ){


		Aii = A + j      + j*lda;
		Qii = Q + j      + j*ldq;
		Tki = T + (j%mt) + j*ldt;

//		for (i=0;i<vb;i++) tau[i]=Tki[(i%mt)+i*ldt];
//		info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml, vb, Aii, lda, Qii, ldq ); 
//		info = LAPACKE_dorgqr_work( LAPACK_COL_MAJOR, ml, vb, vb, Qii, ldq, tau, work, lwork1 );

//		info = lila_dormqrbz_w03 ( ml, vb, k, i, j, mt, A, lda, Q, ldq, T, ldt, work, lwork1 );	

		k  += vb;
		j  += vb;
		l  += vb;
		ml -= vb;

		if( j + mt >= i + n ) vb = n - ( j - i ); else vb = mt;

	}

	return 0;
}

