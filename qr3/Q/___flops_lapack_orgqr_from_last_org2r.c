#include "qr3.h"

long int flops_lapack_orgqr_from_last_org2r( int int_m, int int_n, int int_k, int int_b ){

	long int m, k, b, flops;
	long int kb;

	flops = (( long int ) 0 );

	m = ( long int ) int_m;
	k = ( long int ) int_k;
	b = ( long int ) int_b;

	if( k%b == 0 ) kb = k/b - 1; else kb = k/b;

	flops = kb * (  
		+ ( (long int) 3 ) * m * b * b 
		+ ( (long int) 3 ) * (m-b*kb) * b * b 
		+ ( (long int) 1 ) * b * b * b 
		+ ( (long int) 3 ) * m * b
		+ ( (long int) 3 ) * (m-b*kb) * b
		- ( (long int) 1 ) * b
		- ( (long int) 3 )
		) / ( (long int) 3 );

//	flops = kb * (  
//		+ ( (long int) 6 ) * m  * b * b 
//		- ( (long int) 3 ) * kb * b * b * b
//		+ ( (long int) 1 ) * b * b * b 
//		+ ( (long int) 6 ) * m * b
//		- ( (long int) 3 ) * kb * b  * b
//		- ( (long int) 1 ) * b
//		- ( (long int) 3 )
//		) / ( (long int) 3 );


	return flops;

}
