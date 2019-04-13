#include "qr3.h"

unsigned long int flops_ULTinU_check( int n ){

	unsigned long int flops;

	int n1, n2; 

	flops = (( unsigned long int ) 0 );

	if ( n <= 1 ){

	} else { 

		n1 = n/2;

		n2 = n-n1;

		flops += flops_trmm( 'R', n1, n2 );

		flops += flops_ApUBTinA ( n1, n2 );

		flops += flops_ULTinU_check( n1 );

		flops += flops_ULTinU_check( n2 );
		
	}

	return flops;

}
