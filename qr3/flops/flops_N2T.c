#include "flops.h"

long int flops_N2T( int int_n ){

	long int n, flops;

	n = ( long int ) int_n;

	flops = ( n * n * n - n ) / 3;

	return flops;

}
