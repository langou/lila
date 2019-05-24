#include "flops.h"

long int flops_geqr3_ISW_UT( int int_m, int int_n ){

	long int m, n, flops;
	
	m = (( long int ) int_m );
	n = (( long int ) int_n );

	flops = (( long int ) 0 );

	// flops of geqr3_ISW
	flops = (
	          98 * m * n * n 
	       -  25 * n * n * n 
	       +  28* m 
	       + 217 * n 
	       + 18 
	) / 42

	// This is the full recursion saving term by using UT in the qr3 framework
	// But for ISW we only save some of these flops -- during the connect when we call to the right
	- (n*n*n-n)/6;
	
	// How to find the recursion relation

	return flops;

}
