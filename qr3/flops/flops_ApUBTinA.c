#include "flops.h"

long int flops_ApUBTinA( int int_m, int int_n ){

	long int m, n, flops;

	m = ( long int ) int_m;

	n = ( long int ) int_n;

	flops = m * m * n + m * n ;

	return flops;

}
