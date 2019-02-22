#include "lila.h"

int lila_dgeqrf_q03_3( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork ){

	double *Ajj, *Qjj, *Tkj, *tau;
	int vb, info, j, l, ml, k, lwork1;
//
//	This code block is working from the bottom right of the updated matrix and working by mt - or trying
//

	tau  = work;
	work = work + n;
	lwork1 = lwork - n;
 
	vb = mt - ( (n+i)%mt ); if ( vb > n ) vb = n;
	j  = n + i;
	ml = m - n - i;
	k  = n;	

	Ajj = A + j      + j*lda;
	Qjj = Q + j      + j*ldq;
	Tkj = T + (j%mt) + j*ldt;

	for (i=0;i<vb;i++) tau[i]=Tkj[(i%mt)+i*ldt];
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml, vb, Ajj, lda, Qjj, ldq ); 
	info = LAPACKE_dorgqr_work( LAPACK_COL_MAJOR, ml, vb, vb, Qjj, ldq, tau, work, lwork1 );
	info = lila_dormqrbz_w03( m, vb, k, i, j, mt, A, lda, Q, ldq, T, ldt, work, lwork1 );	

	while( vb != i ){


		if( j - mt <= i ) vb = i; else vb = mt;
		k  -= vb;
		j  -= vb;
		ml += vb;

		Ajj = A + j      + j*lda;
		Qjj = Q + j      + j*ldq;
		Tkj = T + (j%mt) + j*ldt;

		for (i=0;i<vb;i++) tau[i]=Tkj[(i%mt)+i*ldt];
		info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml, vb, Ajj, lda, Qjj, ldq ); 
		info = LAPACKE_dorgqr_work( LAPACK_COL_MAJOR, ml, vb, vb, Qjj, ldq, tau, work, lwork1 );
		info = lila_dormqrbz_w03( m, vb, k, i, j, mt, A, lda, Q, ldq, T, ldt, work, lwork1 );	


	}

/*
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
*/
	return 0;
}

