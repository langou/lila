#include "qr3.h"

long int flops_geqr3_wob_check( int m, int n ){

//	this checks that, indeed, without T-blocking the recursive QR factorization has the same cost as QR2
//	this routine should return the same value as flops_lapack_geqr2 (and it does)
//	the line commented corresponds to the extra flops due to the use of T 
//	if one choses to not use T (not create and not use in LARFB, but use a loop around LARF) 
//	then one would get this cost

	long int flops;

	int n1, n2;

	flops = (( long int ) 0 );

	if ( n == 1){

		flops += flops_lapack_larfg( m );

	}
	else {

		n1 = n/2;
		n2 = n-n1;

		flops += flops_geqr3_wob_check( m, n1 );

		flops += n1 * n1 * n2 ;            // 1: TRMM
		flops += 2 * n1 * n2 * ( m-n1 );   // 2: GEMM
//		flops += (n1-1) * n1 * n2 ;        // 3: TRMM (extra from a bunch of LARF)
		flops += n1 * n2 ;                 // 3: TRMM
		flops += 2 * n1 * n2 * ( m-n1 );   // 4: GEMM
		flops += n1 * n1 * n2;             // 5: TRMM

		flops += flops_geqr3_wob_check( m-n1, n2 );

//		flops += n1 * n2 * n2 ;            // 1: TRMM
//		flops += n1 * n2 * (m-n) ;         // 2: GEMM
//		flops += n1 * n1 * n2 ;            // 3: TRMM
//		flops += n1 * n2 * n2 ;            // 4: TRMM

	}

	return flops;

}
