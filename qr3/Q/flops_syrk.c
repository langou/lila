#include "qr3.h"

long int flops_syrk( int int_n, int int_k ){

	long int n, k, flops;

	n = ( long int ) int_n;
	k = ( long int ) int_k;

	flops = k * n * n;

	return flops;

}
