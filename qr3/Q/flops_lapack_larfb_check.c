#include "qr3.h"

long int flops_lapack_larfb_check( int m, int n, int k ){

	long int flops;

	flops = (( long int ) 0 );

	flops += flops_trmm( 'L', k, n );

	flops += flops_gemm( k, n, m-k );

	flops += flops_trmm( 'L', k, n );

	flops += flops_gemm( m-k, n, k );

	flops += flops_trmm( 'L', k, n );

//	flops += (( long int ) n ) * (( long int ) k );

	return flops;
	
}
