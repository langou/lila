#include "qr3.h"

long int flops_lapack_orgqr( int int_m, int int_n, int int_k, int int_b ){

	long int m, n, k, b, flops;

	flops = (( long int ) 0 );

	m = ( long int ) int_m;
	n = ( long int ) int_n;
	k = ( long int ) int_k;
	b = ( long int ) int_b;

	flops = ( (( long int ) 12 ) * m * n * k
		- (( long int )  6 ) * m * k * k 
		- (( long int )  6 ) * n * k * k 
		+ (( long int )  4 ) * k * k * k
		+ (( long int )  6 ) * m * k 
		- (( long int )  3 ) * k * k 
		- (( long int )  1 ) * k 
		- (( long int )  3 ) ) 
		/ (( long int )  3 );

	return flops;

}
