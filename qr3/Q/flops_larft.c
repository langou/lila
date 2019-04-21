#include "qr3.h"

long int flops_larft( int int_m, int int_k ){

	long int m, k, flops;

	m = ( long int ) int_m;
	k = ( long int ) int_k;

	flops = ((( long int ) 3) * m * k * k - k * k * k - (( long int ) 3) * m + (( long int ) 3) ) / (( long int ) 3) ;

	return flops;

}
