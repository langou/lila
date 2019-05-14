#include "flops.h"

long int flops_lapack_orgqr_bef_check( int m, int n, int k, int nb ){

	int k0, m1, n1, n2, ib;

	long int flops;

	flops = (( long int ) 0 );
	
	ib = k%nb; if (ib==0) ib=nb;

	k0 = k-ib;

	m1 = m-k0;
	n1 = n-k0;

	while( k0 > 0 ){

		ib = nb; if( k0 - ib < 0 ) ib = k0;

		n2 = n1;
		m1 += ib;		
		n1 += ib;
		k0 -= ib;		
	
		flops += flops_larft( m1, ib );

		flops += (ib-1) * ib * n2 ;        // 3: TRMM (extra from a bunch of LARF)

	}

	return flops;

}
