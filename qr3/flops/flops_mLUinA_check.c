#include "flops.h"

long int flops_mLUinA_check( int n ){

	long int flops;

	int n1, n2;

	flops = (( long int ) 0 );

	if ( n <= 1 ){

	} else {

		n1 = n/2; 

		n2 = n - n1;

		flops += flops_mLUinA( n2 );

		flops += flops_gemm( n2, n2, n1 );

		flops += flops_trmm( 'L', n1, n2 );

		flops += flops_trmm( 'R', n2, n1 );

		flops += flops_mLUinA( n1 );

	}

	return flops;

}
