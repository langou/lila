#include "flops.h"

long int flops_geqr3_ISW_constructT( int int_m, int int_n ){

//	this formula is only exact for n being a power of two
//	this formula is a good approximation otherwise
//	flops = 1/3 mn^2 - mn + 2/3 m - 2/21 n^3  + 1/2 n^2 - 5/6 n + 3/7  ;

	long int flops, m, n;

	m = (( long int ) int_m );
	n = (( long int ) int_n );

	flops = (14 * m*n*n - 42 * m*n + 28* m - 4*n*n*n  +21 *n*n -35* n + 18 )/42  ;
	
	return flops;

}
