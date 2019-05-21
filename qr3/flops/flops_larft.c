#include "flops.h"

long int flops_larft( int int_m, int int_n ){

	long int m, n, flops;

	m = ( long int ) int_m;
	n = ( long int ) int_n;

	flops = ( n * ( n - 1 ) * ( 6 * m - ( 2 * n - 1 ) )  ) / 6;

	return flops;

}
