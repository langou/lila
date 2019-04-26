#include "qr3.h"

long int flops_lapack_larfb( int int_m, int int_n, int int_k ){

	long int m, n, k, flops;

	m = ( long int ) int_m;
	n = ( long int ) int_n;
	k = ( long int ) int_k;

	flops = (( long int ) 0 );

	flops += 4 * m * n * k - n * k * k + 2 * n * k;

	return flops;

}

//	the flops count is "tweaked" so that the cost of LARFB( m, n, k )
//	is the same as the flops count of 
//		for i = 1:k,
//	  		 LARF( m-i+1, n )
//		end
//	and the same as of 
//		for i = 1:k,
//			LARFB( m-i+1, 1 )
//		end

//	instead of doing:
//	flops += k * k * n ;            // 1: TRMM
//	flops += 2 * k * n * ( m-k );   // 2: GEMM
//	flops += k * k * n ;            // 3: TRMM
//	flops += 2 * k * n * ( m-k );   // 4: GEMM
//	flops += k * k * n ;            // 5: TRMM
//	flops += n * k;                 // 6: ADD

//	we do:
//	flops += k * k * n + n * k ;    // 1: TRMM
//	flops += 2 * k * n * ( m-k );   // 2: GEMM
//	flops += k * k * n ;            // 3: TRMM
//	flops += 2 * k * n * ( m-k );   // 4: GEMM
//	flops += k * k * n + n * k ;    // 5: TRMM
//	flops += 0;                     // 6: ADD
//
//	why
//
//	step 1 and 4, we add n * k. We use the FLOPS count of TRMM as n * k * k
//	we assume a "smart" inner product that does 2*n-1 operations, so we say
//	n multiplications and (n-1) additions. However when we count FLOPS in other
//	operation, we do not do this trick. To be consistent with LARF, we count
//	these two TRMM as n * k * k + n * k (which means that we use 2*n FLOPS inner products)
//
//	step 3, for step 3 we count TRMM n * k * k, this is because we really want this TRMM
//	to only count as n in the case k=1. (n * k * k + n * k would give 2*n).
//	LARF just does a scaling in the case k=1
//
//	the fact that we do not treat ( TRMM 1 and 5 ) and ( TRMM 3 ) not the same makes sense
//	( TRMM 3 ) is a multiplication
//	( TRMM 1 and 5 ) is part of an inner product that is broken in two pieces,
//	so we need to do the multiplication but we still have to an addition to accumulate
//	the pieces
//
//	Note that we do not take into account whether matrices are unit triangular or triangular.
//	This is because LARF does not take into account whether there is a 1 at the top of v_i
//
//	The step 6 is zero because we consider that we do the + at the same time as TRMM 5.
