#include "qr3.h"

long int flops_lapack_org2r( int int_m, int int_n, int int_k ){

	long int m, n, k, flops;

	m = ( long int ) int_m;
	n = ( long int ) int_n;
	k = ( long int ) int_k;

	flops = ( (( long int ) 12 ) * m * n * k
		- (( long int )  6 ) * m * k * k 
		- (( long int )  6 ) * n * k * k 
		+ (( long int )  4 ) * k * k * k
		+ (( long int )  9 ) * n * k 
		- (( long int )  3 ) * m * k 
		- (( long int )  3 ) * k * k
		- (( long int )  4 ) * k ) 
		/ (( long int )  3 );

	return flops;

}
