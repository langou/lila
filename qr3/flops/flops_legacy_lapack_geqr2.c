#include "flops.h"

long int flops_legacy_lapack_geqr2( int int_m, int int_n ){

// this returns the number of FLOPS as per LAPACK Working Note #41
// the formula in LAPACK Working Note #41 is
// (2)mn^2 - (2/3)(n^3) + (2)(mn) + (17/3)(n)

	long int m, n, flops;

	m = ( long int ) int_m;
	n = ( long int ) int_n;

	flops = ( (( long int ) 6) * m * n * n 
	- (( long int ) 2) * n * n * n 
	+ (( long int ) 6) * m * n 
	+ (( long int ) 17) * n )
	/ (( long int ) 3);

	return flops;

}
