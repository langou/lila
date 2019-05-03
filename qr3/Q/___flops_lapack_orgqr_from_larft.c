#include "qr3.h"

long int flops_lapack_orgqr_from_larft( int int_m, int int_n, int int_k, int int_b ){

	long int m, k, b, flops;
	long int kb;

	flops = (( long int ) 0 );

	m = ( long int ) int_m;
	k = ( long int ) int_k;
	b = ( long int ) int_b;

	if( k%b == 0 ) kb = k/b - 1; else kb = k/b;

	flops = 
	  (( long int ) 6) * ( m * kb * b * b ) 
	- (( long int ) 6) * ( kb * kb * b * b * b ) 
	+ (( long int ) 6) * ( kb * b * b * b * ( kb + 1 ) / 2 )
	- (( long int ) 6) * ( m * kb * b ) 
	- (( long int ) 6) * ( kb * kb * b * b ) 
	+ (( long int ) 6) * ( kb * b * b * ( kb + 1 ) / 2 ) 
	- (( long int ) 2) * kb * b * b * b
	+ (( long int ) 2) * kb * b * b
	+ kb * b * b
	- kb * b;
	//printf(" check ==== %ld \n",flops%6);
	flops = flops / (( long int ) 6);

	return flops;

}





