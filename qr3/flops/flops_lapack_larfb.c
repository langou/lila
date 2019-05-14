#include "flops.h"

long int flops_lapack_larfb( int int_m, int int_n, int int_k ){

	long int m, n, k, flops;

	m = ( long int ) int_m;
	n = ( long int ) int_n;
	k = ( long int ) int_k;

	flops = (( long int ) 0 );

	flops += k * k * n ;            // 1: TRMM
	flops += 2 * k * n * ( m-k );   // 2: GEMM
	flops += (k-1) * k * n ;        // 3: TRMM (extra from a bunch of LARF)
	flops += k * n ;                // 3: TRMM
	flops += 2 * k * n * ( m-k );   // 4: GEMM
	flops += k * k * n;             // 5: TRMM

	return flops;

}
