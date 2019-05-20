#include "flops.h"

long int flops_geqr3_ISW( int m, int n ){

	long int flops;

	int n1, n2;

	flops = (( long int ) 0 );

	if ( n == 1){

		flops += flops_lapack_larfg( m );

	} else {

		n1 = n/2;
		n2 = n-n1;

		flops += flops_geqr3_check( m, n1 );

		flops += n1 * n1 * n2;                           // trmm
		flops += (( long int ) 2 ) * n1 * n2 * ( m-n1 ); // gemm
		flops += ( n1-1 ) * n1 * n2;                     // trmm
		flops += n1 * n2;                                // trmm
		flops += (( long int ) 2 ) * n1 * n2 * ( m-n1 ); // gemm
		flops += n1 * n1 * n2;                           // trmm

		flops += flops_geqr3_ISW( m-n1, n2 );
//		flops += flops_geqr3_check( m-n1, n2 );

	}

	return flops;

}
