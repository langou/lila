#include "flops.h"

long int flops_geqr3_ISW_UT( int int_m, int int_n ){

	long int m, n, flops;
	
	m = (( long int ) int_m );
	n = (( long int ) int_n );

//	flops = ( 98 * m * n * n + 28 * m - 27*n*n*n + 231 * n + 6 )/42  ;

	flops = flops_geqr3_ISW( int_m, int_n );
	flops -= ( n*n*n - 7 * n + 6 ) /21;

	return flops;

}
