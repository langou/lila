#include "flops.h"

long int flops_larft( int int_m, int int_k ){

	long int m, k, flops;

	m = ( long int ) int_m;
	k = ( long int ) int_k;

	flops = ( k * ( k - (( long int ) 1) ) * ( (( long int ) 6) * m - ( (( long int ) 2) * k - (( long int ) 1) ) )  ) / (( long int ) 6);

	return flops;

}
