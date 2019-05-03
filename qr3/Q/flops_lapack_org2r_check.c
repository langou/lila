#include "qr3.h"

long int flops_lapack_org2r_check( int m, int n, int k ){

	int k0, m1, n1, n2;

	long int flops;

	flops = (( long int ) 0 );
	
	k0 = k-1;

	m1 = m-k0;
	n1 = n-k0;

//	dorg2r_( &m1, &n1, &ib, A11, &lda, tau1, work, &lwork );
//	flops += flops_lapack_org2r( m1, n1, ib );

	flops += flops_lapack_larfb( m1, n1, 1 );
	flops += 4 * m1 - 3;

	while( k0 > 0 ){

		n2 = n1;
		m1 ++;		
		n1 ++;
		k0 --;		
	
//		LAPACKE_dlarft_work( LAPACK_COL_MAJOR, 'F', 'C', m1, ib, A11, lda, tau1, work, ib);
//		flops += flops_larft( m1, ib );

//		our_dlarfb_lnfc( m1, n2, ib, A11, lda, work, ib, A12, lda, work+ib*ib );
		flops += flops_lapack_larfb( m1, n2, 1 );

//		dorg2r_( &m1, &ib, &ib, A11, &lda, tau1, work, &lwork );
//		flops += flops_lapack_org2r( m1, 1, 1 );
		flops += flops_lapack_larfb( m1, 1, 1 );
		flops += 4 * m1 - 3;

	}	

	return flops;

}
