#include "qr3.h"

long int flops_N2T( int int_n ){

	long int n, flops;

	n = ( long int ) int_n;

	flops = ( n * n * n - n ) / (( long int ) 3 );

	return flops;

}
