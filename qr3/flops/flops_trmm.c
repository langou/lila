#include "flops.h"

long int flops_trmm( char S, int int_m, int int_n ){

	long int m, n, flops;

	m = ( long int ) int_m;
	n = ( long int ) int_n;

	if( S == 'L' ){
		flops = n * m * m;
	}

	if( S == 'R' ){
		flops = m * n * n;
	}

	return flops;

}
