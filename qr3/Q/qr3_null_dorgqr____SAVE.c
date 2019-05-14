#include "qr3.h"

int qr3_null_dorgqr( int m, int n, int k, int nb, double *A, int lda, double *tau, double *work, int lwork ){

	double *A11, *Q11, *tau1;
	int k0, k1, m1, n1, n2, ib;
	
	ib = k%nb; if (ib==0) ib=nb;

	k0 = k-ib;
	k1 = n-ib;

	m1 = m-k0;
	n1 = n-k0;
	//m1 = m-k;
	//n1 = n-k;

	A11 = A+k0*(1+lda);
	Q11 = A+(k-ib)+k*lda;
	//A11 = A;
	//Q11 = A+k+k*lda;
	tau1 = tau+k0;

	dorgqr_after( m1, n, ib, A11, lda, A11, lda, Q11, lda );
	//dorgqr_after( m, n, k, A11, lda, A11, lda, Q11, lda );

	while( k0 > 0 ){

		ib = nb; if( k0 - ib < 0 ) ib = k0;

		A11  -= ib*(1+lda);
		Q11  -= ib;
		tau1 -= ib;

		n2 = n1;
		m1 += ib;		
		n1 += ib;
		k0 -= ib;		

		dorgqr_after( m1, n, ib, A11, lda, A11, lda, Q11, lda );
	
		//LAPACKE_dlarft_work( LAPACK_COL_MAJOR, 'F', 'C', m1, ib, A11, lda, tau1, work, ib);
		//our_dlarfb_lnfc( m1, n-k, ib, A11, lda, work, ib, Q11, lda, work+ib*ib );

//		our_dlarfb_lnfc( m1, n-k, ib, A11, lda, A11, lda, Q11, lda, work );

	}	

	return 0;

}
