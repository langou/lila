#include "qr3.h"

int our_dorgqr( int m, int n, int k, int nb, double *A, int lda, double *tau, double *work, int lwork, int info ){

	double *Akk, *A0k, *Aik, *tauk;
	int kk, ml, nl, ib, i, j, ldwork;
	
	ldwork = n;

	kk    = n;

	ml  = m-kk;
	nl  = n-kk;

	A0k  = A+kk*lda;
	Akk  = A+kk+kk*lda;
	tauk = tau+kk;

	ib = nb; if( kk - ib < 0 ) ib = n - nl;

	A0k  -= ib*lda;	
	Akk  -= ib*(1+lda);
	tauk -= ib;

	ml += ib;		
	nl += ib;		
	kk -= ib;

	for( i = 0; i < kk; i++){ for( j = 0; j < ib; j++ ){ A0k[i+j*lda] = (+0.0e00); } }

	dorg2r_( &ml, &ib, &ib, Akk, &lda, tauk, work, &lwork );

	while( kk > 0 ){

		ib = nb; if( kk - ib < 0 ) ib = n - nl;

		A0k  -= ib*lda;
		Akk  -= ib*(1+lda);
		Aik   = Akk+ib*lda;
		tauk -= ib;

		ml += ib;		
		nl += ib;		
		kk -= ib;
	
//		info = LAPACKE_dlarft_work( LAPACK_COL_MAJOR, 'F', 'C', ml, ib, Akk, lda, tauk, work, ldwork);
		info = LAPACKE_dlarft_work( LAPACK_COL_MAJOR, 'F', 'C', ml, ib, Akk, lda, tauk, work, ib);

//		info = LAPACKE_dlarfb_work( LAPACK_COL_MAJOR, 'L', 'N', 'F', 'C', ml, nl-ib, ib, Akk, lda, work, ldwork, Akk+ib*lda, lda, work+ib, ldwork);
		info = our_dlarfb_lnfc( ml, nl-ib, ib, Akk, lda, work, ib, Akk+ib*lda, lda, work+ib*ib );

		dorg2r_( &ml, &ib, &ib, Akk, &lda, tauk, work, &lwork );

		for( i = 0; i < kk; i++){ for( j = 0; j < ib; j++ ){ A0k[i+j*lda] = (+0.0e00); } }

	}	

	return 0;

}
