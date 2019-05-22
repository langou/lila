#include "flops.h"

//	flops = mn^2 - 1/3 n^3 - mn + 1/2 n^2 - 1/6 n

long int flops_larft( int int_m, int int_n ){

	long int m, n, flops;

	m = ( long int ) int_m;
	n = ( long int ) int_n;

//	flops = ( n * ( n - 1 ) * ( 6 * m - ( 2 * n - 1 ) )  ) / 6;

	flops = ( 6 * m * n * n
	        - 2 * n * n * n
	        - 6 * m * n 
	        + 3 * n * n 
	        - 1 * n ) / 6;

	return flops;

}
