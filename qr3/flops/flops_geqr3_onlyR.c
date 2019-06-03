#include "flops.h"

long int flops_geqr3_onlyR( int int_n ){

	long int flops;
	long int n1, n2;
	long int n;

	n = ( long int ) int_n;
	flops = (( long int ) 0 );

//	if ( n == 1){
//
//	} else {
//
//	n1 = n/2;
//	n2 = n-n1;
//
//	flops += flops_geqr3_onlyR_check( n1 );
//	flops += flops_geqr3_onlyR_check( n2 );
//
//	flops += flops_trmm( 'L', n1, n2 );
//
//	}

//	This is flops for USE_T
//	( n * n * n - 3 * n * n + 2 * n ) / 6

	return flops;

}
