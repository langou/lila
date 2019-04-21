#include "qr3.h"

long int flops_gemm( int int_m, int int_n, int int_k ){

	long int m, n, k, flops;

	m = ( long int ) int_m;
	n = ( long int ) int_n;
	k = ( long int ) int_k;

	flops = (( long int ) 2) * m * n * k;

	return flops;

}


