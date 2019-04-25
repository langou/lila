#include "qr3.h"

long int flops_lapack_geqr2( int int_m, int int_n ){

	long int m, n, flops;

	m = ( long int ) int_m;
	n = ( long int ) int_n;

	// formula for GEQRF based off of working note 41 - no nb included   ---- I multiplied 3 through so no divisions are done
	flops = ( (( long int ) 6) * m * n * n 
	- (( long int ) 2) * n * n * n 
	+ (( long int ) 3) * m * n 
	+ (( long int ) 3) * n * n 
	+ (( long int ) 14) * n )
	/ (( long int ) 3);

	return flops;

}
