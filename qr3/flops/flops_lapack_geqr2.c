#include "flops.h"

long int flops_lapack_geqr2( int int_m, int int_n ){

	long int m, n, flops;

	m = ( long int ) int_m;
	n = ( long int ) int_n;

	flops = ( 6 * m * n * n 
	       -  2 * n * n * n 
	       +  3 * m * n 
	       + 17 * n ) /  3;

	return flops;

}
