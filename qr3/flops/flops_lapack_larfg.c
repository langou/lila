#include "flops.h"

long int flops_lapack_larfg( int int_m ){

	long int m, flops;

	m = ( long int ) int_m;

	flops = (( long int ) 3) * m + 5;

	return flops;

}

