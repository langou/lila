#include "flops.h"

long int flops_mLUinA( int int_n ){

	long int n, flops;

	n = ( long int ) int_n;

	flops = ( (( long int ) 2 ) * ( n * n * n - n ) ) / (( long int ) 3 );

	return flops;

}
