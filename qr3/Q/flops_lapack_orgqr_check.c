#include "qr3.h"

long int flops_lapack_orgqr_check( int m, int n, int k, int nb ){

	long int flops;

	int k0, m1, n1, n2, ib;

	flops = (( long int ) 0 );
	
	ib = k%nb; if (ib==0) ib=nb;

	k0 = k-ib;

	m1 = m-k0;
	n1 = n-k0;

	//flops += flops_org2r( m1, n1, ib );

	int kb;
	kb = 1;

	while( k0 > 0 ){

		ib = nb; if( k0 - ib < 0 ) ib = k0;

		n2 = n1;
		m1 += ib;		
		n1 += ib;
		k0 -= ib;		

//		flops += flops_larft( m1, ib );
	
//		flops += flops_lapack_larfb( m1, n2, ib );

//		flops += flops_org2r( m1, ib, ib );

	}	

	return flops;

}
