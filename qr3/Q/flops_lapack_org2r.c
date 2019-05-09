#include "qr3.h"

long int flops_lapack_org2r( int int_m, int int_n, int int_k ){

// this returns the number of FLOPS for ORG2R as per our computation
// (4)mnk - (2)(m+n)(k^2) + (4/3)(k^3) - (mk) + (nk) - (1/3)(k)

// the number of FLOPS for ORG2R as per LAPACK Working Note #41 is
// (4)mnk - (2)(m+n)(k^2) + (4/3)(k^3) - (mk) + (3)(nk) - (1)(k^2) - (4/3)(k)
// difference with ~ ............................~ ........~ .........~......
// our formula:
// (4)mnk - (2)(m+n)(k^2) + (4/3)(k^3) - (mk) + (1)(nk) - (0)(k^2) - (1/3)(k) 

	long int m, n, k, flops;

	m = ( long int ) int_m;
	n = ( long int ) int_n;
	k = ( long int ) int_k;
	
	flops = ( (( long int ) 12 ) * m * n * k
		- (( long int )  6 ) * m * k * k
		- (( long int )  6 ) * n * k * k
		+ (( long int )  4 ) * k * k * k
		- (( long int )  3 ) * m * k
		+ (( long int )  3 ) * n * k
		- (( long int )  1 ) * k )
		/ (( long int )  3 );

	return flops;

}
