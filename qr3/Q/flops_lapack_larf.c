#include "qr3.h"

long int flops_lapack_larf( int int_m, int int_n ){

	long int m, n, flops;

	m = ( long int ) int_m;
	n = ( long int ) int_n;

	flops = (( long int ) 4) * m * n - n;

	return flops;

}
