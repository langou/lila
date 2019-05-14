#include "flops.h"

long int flops_V2N( int int_n ){

	long int n, flops;

	n = ( long int ) int_n;

	flops = ( n * n * n - n ) / (( long int ) 3 );

	return flops;

}
