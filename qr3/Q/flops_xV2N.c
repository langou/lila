#include "qr3.h"

unsigned long int flops_xV2N( int n ){

	int n1, n2;
	unsigned long int flops;

	flops = 0;

	if ( n <= 1 ) {

		return 0;


	} else {

		n1 = n / 2;
		n2 = n - n1;

		flops += flops_xV2N( n1 );
		flops += flops_syrk( n1, n2 ); // 1/6 flops
		flops += flops_trmm( n1, n2, 'R' ); // 1/3 flops
		flops += flops_xV2N( n2 );

		return flops;

	}

}


