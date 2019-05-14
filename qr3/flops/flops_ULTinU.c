#include "flops.h"

long int flops_ULTinU( int int_n ){

	long int n, flops;

	n = ( long int ) int_n;

	flops = ( (( long int ) 2 ) * n * n * n + (( long int ) 3 ) * n * n - (( long int ) 5 ) * n ) / (( long int ) 6 );

	return flops;

}
