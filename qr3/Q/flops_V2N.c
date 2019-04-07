#include "qr3.h"

unsigned long int flops_V2N( int n ){

	unsigned long int u_n, flops;

	u_n = ( unsigned long int ) n;

	flops = ( u_n * u_n * u_n - u_n ) / (( unsigned long int ) 3 );

	return flops;

}

/*

	int n1, n2;

	flops = 0;

	if ( n <= 1 ) {

		return 0;


	} else {

		n1 = n / 2;
		n2 = n - n1;

		flops += flops_V2N( n1 );
		flops += flops_syrk( n1, n2 );
		flops += flops_trmm( n1, n2, 'R' );
		flops += flops_V2N( n2 );

		return flops;

	}

}

*/
