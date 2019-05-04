#include "qr3.h"

long int flops_lapack_orgqr_check( int m, int n, int k, int nb ){

	int k0, m1, n1, n2, ib;

	long int flops;

	flops = (( long int ) 0 );
	
	ib = k%nb; if (ib==0) ib=nb;

	k0 = k-ib;

	m1 = m-k0;
	n1 = n-k0;

//	dorg2r_( &m1, &n1, &ib, A11, &lda, tau1, work, &lwork );
//	flops += flops_lapack_org2r_check( m1, n1, ib );

	flops += flops_lapack_larfb( m1, n1-ib, ib );
	flops += flops_lapack_org2r_check( m1, ib, ib );

	while( k0 > 0 ){

		ib = nb; if( k0 - ib < 0 ) ib = k0;

		n2 = n1;
		m1 += ib;		
		n1 += ib;
		k0 -= ib;		
	
//		LAPACKE_dlarft_work( LAPACK_COL_MAJOR, 'F', 'C', m1, ib, A11, lda, tau1, work, ib);
//		flops += flops_larft( m1, ib );

//		our_dlarfb_lnfc( m1, n2, ib, A11, lda, work, ib, A12, lda, work+ib*ib );
		flops += flops_lapack_larfb( m1, n2, ib );

//		dorg2r_( &m1, &ib, &ib, A11, &lda, tau1, work, &lwork );
		flops += flops_lapack_org2r_check( m1, ib, ib );

	}

	return flops;

}
