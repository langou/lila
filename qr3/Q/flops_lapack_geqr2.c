#include "qr3.h"

long int flops_lapack_geqr2( int int_m, int int_n ){

	long int m, n, flops;

	m = ( long int ) int_m;
	n = ( long int ) int_n;

	flops = ( (( long int ) 6) * m * n * n 
	- (( long int ) 2) * n * n * n 
	+ (( long int ) 3) * m * n 
	+ (( long int ) 17) * n )
	/ (( long int ) 3);

	return flops;

}
