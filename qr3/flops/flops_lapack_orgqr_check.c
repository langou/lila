#include "flops.h"

long int flops_lapack_orgqr_check( int m, int n, int k, int nb ){

	int k0, m1, n1, n2, ib;

	long int flops;

	flops = (( long int ) 0 );
	
	ib = k%nb; if (ib==0) ib=nb;

	k0 = k-ib;

	m1 = m-k0;
	n1 = n-k0;

//	flops += flops_lapack_org2r_check( m1, n1, ib );
	flops += flops_lapack_org2r( m1, n1, ib );

	while( k0 > 0 ){

		ib = nb; if( k0 - ib < 0 ) ib = k0;

		n2 = n1;
		m1 += ib;		
		n1 += ib;
		k0 -= ib;		
	
		flops += flops_larft( m1, ib );

		flops += flops_lapack_larfb( m1, n2, ib );

//		flops += ib * ib * n2 ;             // 1: TRMM (these flops are saved during a BZ)
//		flops += 2 * ib * n2 * ( m1-ib );   // 2: GEMM
//		flops += (ib-1) * ib * n2 ;         // 3: TRMM (extra from a bunch of LARF)
//		flops += ib * n2 ;                  // 3: TRMM
//		flops += 2 * ib * n2 * ( m1-ib );   // 4: GEMM
//		flops += ib * ib * n2;              // 5: TRMM

//		flops += flops_lapack_org2r_check( m1, ib, ib );
		flops += flops_lapack_org2r( m1, ib, ib );

	}

	return flops;

}
