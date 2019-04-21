#include "qr3.h"

long int flops_V2N_check( int n ){

	long int flops;

	int n1, n2;

	flops = (( long int ) 0 );

	if ( n <= 1 ) {

		return (( long int ) 0 );


	} else {

		n1 = n / 2;
		n2 = n - n1;

		flops += flops_V2N( n1 );

		flops += flops_syrk( n1, n2 );

		flops += flops_trmm( 'R', n1, n2 );

		flops += flops_V2N( n2 );

		return flops;

	}

}

