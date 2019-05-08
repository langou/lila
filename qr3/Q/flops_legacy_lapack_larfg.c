#include "qr3.h"

long int flops_legacy_lapack_larfg( int int_m ){

	long int m, flops;

	m = ( long int ) int_m;

	flops = (( long int ) 4) * m + 5;

	return flops;

}

