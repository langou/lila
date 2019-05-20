#include "flops.h"

long int flops_geqr3_ISW_save( int m, int n ){

	long int flops;

	int n1, n2;

	flops = (( long int ) 0 );

	if ( n == 1){
		//Do nothing
	} else {

		n1 = n/2;
		n2 = n-n1;

		flops += flops_geqr3_ISW_save( m, n1 );

//		flops += n1 * n2 * n2 ;            // 1: TRMM
//		flops += n1 * n2 * (m-n) ;         // 2: GEMM
//		flops += n1 * n1 * n2 ;            // 3: TRMM
//		flops += n1 * n2 * n2 ;            // 4: TRMM

		flops += n1 * n2 * n2;
		flops += n1 * n2 * m;
		
		flops += flops_geqr3_ISW_save( m-n1, n2 );

	}

	return flops;

}
