#include "qr3.h"

unsigned long int flops_lapack_larfb_check( int m, int n, int k ){

	unsigned long int flops;

	flops = (( unsigned long int ) 0 );

	flops += flops_trmm( 'L', k, n );

	flops += flops_gemm( k, n, m-k );

	flops += flops_trmm( 'L', k, n );

	flops += flops_gemm( m-k, n, k );

	flops += flops_trmm( 'L', k, n );

	flops += (( unsigned long int ) n ) * (( unsigned long int ) k );

	return flops;

}
