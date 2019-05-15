#include "qr3.h"

int our_dorgqr_Q2( int m, int n, int k, int nb, double *A, int lda, double *tau, double *work, int lwork ){

	double *A11, *A01, *A12, *tau1;
	int k0, m1, n1, n2, ib, i, j, ldwork;
	
	ldwork = n;

	ib = k%nb; if (ib==0) ib=nb;

	k0 = k-ib;

	m1 = m-k0;
	n1 = n-k0;

	A01 = A+k0*lda;
	A11 = A+k0*(1+lda);
	tau1 = tau+k0;

	our_dorg2r_Q2( m1, n1, ib, A11, lda, tau1, work, lwork );
//	dorg2r_( &m1, &n1, &ib, A11, &lda, tau1, work, &lwork );

	for( i = 0; i < k0; i++){ for( j = 0; j < n1; j++ ){ A01[i+j*lda] = (+0.0e00); } }

	while( k0 > 0 ){

		ib = nb; if( k0 - ib < 0 ) ib = k0;

		A01  -= ib*lda;
		A11  -= ib*(1+lda);
		A12   = A11+ib*lda;
		tau1 -= ib;

		n2 = n1;
		m1 += ib;		
		n1 += ib;
		k0 -= ib;		
	
//		LAPACKE_dlarft_work( LAPACK_COL_MAJOR, 'F', 'C', m1, ib, A11, lda, tau1, work, ldwork);
		LAPACKE_dlarft_work( LAPACK_COL_MAJOR, 'F', 'C', m1, ib, A11, lda, tau1, work, ib);

//		LAPACKE_dlarfb_work( LAPACK_COL_MAJOR, 'L', 'N', 'F', 'C', m1, n2, ib, A11, lda, work, ldwork, A12, lda, work+ib, ldwork);
		our_dlarfb_lnfc( m1, n2, ib, A11, lda, work, ib, A12, lda, work+ib*ib );

		dorg2r_( &m1, &ib, &ib, A11, &lda, tau1, work, &lwork );

		for( i = 0; i < k0; i++){ for( j = 0; j < ib; j++ ){ A01[i+j*lda] = (+0.0e00); } }

	}	

	return 0;

}
