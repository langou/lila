#include "flops.h"

long int flops_geqr3_ISW_check( int m, int n ){

	long int flops;

	int n1, n2;

	flops = (( long int ) 0 );

	if ( n == 1){
		flops += flops_lapack_larfg( m );
	} else {

		n1 = n/2;
		n2 = n-n1;

		flops += flops_geqr3_check( m, n1 );

		flops += flops_trmm( 'L', n1, n2 );
		flops += flops_gemm( n1, n2, m-n1 );
		flops += flops_trmm( 'L', n1, n2 );
		flops += flops_gemm( m-n1, n2, n1 );
		flops += flops_trmm( 'L', n1, n2 );

//		for (j=0;j<n2;j++) for (i=0;i<n1;i++) A[(j+n1)*lda+i] -= W[j*ldw+i];
	//	flops += n1 * n2;

		flops += flops_geqr3_ISW_check( m-n1, n2 );
//		flops += flops_geqr3_check( m-n1, n1 );


	}

	return flops;

}
