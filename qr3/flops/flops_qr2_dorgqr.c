#include "flops.h"

long int flops_qr2_dorgqr( int int_m, int int_n ){

	long int flops;
	long int m, n;

	m = ( long int ) int_m;
	n = ( long int ) int_n;

	flops = (( long int ) 0 );

//      this is the cost for creating Q1
	flops = ( 
	        + 6 * (m - n) * n
	        + 4 * n * n
		- 3 * (m - n)
		- 1 )
		* n
		/ (( long int )  3 );

//      this is the cost for creating Q2 --> Note this is zero since k == n
	flops +=  ( 4 * m - 2 * n + 1 ) * n * ( n - n ) ;

//	+ useT - saveBZ
	flops -= n * ( n - 1 ) / (( long int) 2);


	return flops;

}
