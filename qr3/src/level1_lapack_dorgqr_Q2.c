#include "qr2.h"

int level1_lapack_dorgqr_Q2( int m, int n, int k, int nb, double *A, int lda, double *tau, double *work, int lwork ){

	double *A11, *tau1;
	int k0, m1, n1, ib, i, j, ldwork;
	double *Axx, *A0x;
	int nx;
	
	ldwork = n;

	ib = k%nb; if (ib==0) ib=nb;

	k0 = k-ib;

	m1 = m-k0;
	n1 = n-k0;

	A11 = A+k0*(1+lda);
	tau1 = tau+k0;

	A0x = A+k*lda;
	Axx = A+k*(1+lda)-ib;
	nx = n-k;

	for( i = 0; i < k0; i++){ for( j = 0; j < nx; j++ ){ A0x[i+j*lda] = (+0.0e00); } }

	our_dorg2r_Q2( m1, n1, ib, A11, lda, tau1, work, lwork );

	while( k0 > 0 ){

		ib = nb; if( k0 - ib < 0 ) ib = k0;

		A11 -= ib*(1+lda);
		Axx -= ib;
		tau1 -= ib;

		m1 += ib;		
		k0 -= ib;		
	
		LAPACKE_dlarft_work( LAPACK_COL_MAJOR, 'F', 'C', m1, ib, A11, lda, tau1, work, ib);

		our_dlarfb_lnfc( m1, nx, ib, A11, lda, work, ib, Axx, lda, work+ib*ib );

	}	

	return 0;

}
