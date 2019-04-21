#include "qr3.h"

long int flops_lapack_larfb( int int_m, int int_n, int int_k ){

	long int m, n, k, flops;

	m = ( long int ) int_m;
	n = ( long int ) int_n;
	k = ( long int ) int_k;

	flops = (( long int ) 4) * k * m * n - k * k * n - k * n + (( long int ) 1) ;

	return flops;

}
