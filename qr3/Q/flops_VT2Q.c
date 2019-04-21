#include "qr3.h"

long int flops_VT2Q( int int_m, int int_n ){

	long int m, n, flops;

	m = ( long int ) int_m;
	n = ( long int ) int_n;

	flops = m * n * n + ( n * n - n ) / (( long int ) 2 ) ;

	return flops;

}
