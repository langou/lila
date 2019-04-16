#include "qr3.h"

unsigned long int flops_lapack_orgqr_check( int m, int n, int k, int nb ){

	unsigned long int flops;

	int kk, ml, nl, ib;

	flops = (( unsigned long int ) 0 );

	kk    = n;

	ml  = m-kk;
	nl  = n-kk;

	ib = nb; if( kk - ib < 0 ) ib = n - nl;

	ml += ib;		
	nl += ib;		
	kk -= ib;

	flops += flops_org2r( ml, ib, ib );

	while( kk > 0 ){

		ib = nb; if( kk - ib < 0 ) ib = n - nl;

		ml += ib;		
		nl += ib;		
		kk -= ib;
	
		//printf("larft = %lu\n",flops_larft( ml, ib ) - (ml-1));
		flops += flops_larft( ml, ib ) - (ml-1);
		//flops += flops_larft( ml, ib );

		flops += flops_lapack_larfb( ml, nl-ib, ib );

		flops += flops_org2r( ml, ib, ib );

	}	

	return flops;

}
