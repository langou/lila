#include "flops.h"

long int flops_geqr3_UT( int int_m, int int_n ){

	long int m, n, flops;

	m = ( long int ) int_m;
	n = ( long int ) int_n;

//	this will work for n  non-powers of two
//	flops = flops_geqr3_check( int_m, int_n ) - ( n * n * n - n ) / 3;

//	this will only work for powers of two
//	flops = flops_geqr3( int_m, int_n ) - ( n * n * n - n ) / 3;

//	flops = ( 18 * m * n * n -  7 * n * n * n  + 37 * n  ) / 6;

	flops = 3 * m * n * n - n * n * n +6*n - ( n * n * n -  n ) / 6 ;

	return flops;

}
