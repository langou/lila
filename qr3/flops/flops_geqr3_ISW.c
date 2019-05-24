#include "flops.h"

long int flops_geqr3_ISW( int int_m, int int_n ){

	long int m, n, flops;
	
	m = (( long int ) int_m );
	n = (( long int ) int_n );

	flops = (( long int ) 0 );

	flops = (
	          98 * m * n * n 
	       -  25 * n * n * n 
	       +  28* m 
	       + 217 * n 
	       + 18 
	) / 42;

	return flops;

}
