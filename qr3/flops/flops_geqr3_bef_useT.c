#include "flops.h"

long int flops_geqr3_bef_useT( int int_n ){

//	this formula is only exact for n being a power of two
//	this formula is a good approximation otherwise

	long int flops, n;

	flops = (( long int ) 0 );

	n = (( long int ) int_n );

//	flops +=  ( n * ( n - 1 ) * ( n - 2 ) ) / 6;

	flops += ( n * n * n - 3 * n * n + 2 * n ) / 6;
	
	return flops;

}
