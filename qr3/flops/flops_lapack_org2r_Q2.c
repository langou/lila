#include "flops.h"

long int flops_lapack_org2r_Q2( int int_m, int int_n, int int_k ){

// this returns the number of FLOPS for ORG2R as per our computation
// (4)mnk - (2)(m+n)(k^2) + (4/3)(k^3) - (mk) + (nk) - (1/3)(k)

// the number of FLOPS for ORG2R as per LAPACK Working Note #41 is
// (4)mnk - (2)(m+n)(k^2) + (4/3)(k^3) - (mk) + (3)(nk) - (1)(k^2) - (4/3)(k)
// difference with ~ ............................~ ........~ .........~......
// our formula:
// (4)mnk - (2)(m+n)(k^2) + (4/3)(k^3) - (mk) + (1)(nk) - (0)(k^2) - (1/3)(k) 

	long int m, n, k, flops;

	m = ( long int ) int_m;
	n = ( long int ) int_n;
	k = ( long int ) int_k;
	
//	flops = ( (( long int ) 12 ) * m * n * k
//		- (( long int )  6 ) * m * k * k
//		- (( long int )  6 ) * n * k * k
//		+ (( long int )  4 ) * k * k * k
//		- (( long int )  3 ) * m * k
//		+ (( long int )  3 ) * n * k
//		- (( long int )  1 ) * k )
//		/ (( long int )  3 );

//	flops = ( (( long int ) 12 ) * m * n * k
//		- (( long int )  6 ) * m * k * k
//		- (( long int )  6 ) * n * k * k
//		+ (( long int )  4 ) * k * k * k
//		- (( long int ) 12 ) * m * n
//		+ (( long int )  9 ) * m * k
//		+ (( long int ) 15 ) * n * k
//		- (( long int ) 12 ) * k * k
//		- (( long int )  9 ) * n
//		+ (( long int )  8 ) * k )
//		/ (( long int )  3 );
//	flops += 4 * (m-k+1)  * (n-k) - (n-k) ;


//      this is to create Q1
//	flops = ( 
//	        + 6 * (m-k) * k * k
//	        + 4 * k * k * k
//		- 3 * (m-k) * k
//		- 1 * k
//		)
//		/ (( long int )  3 );
//
//      this is the ``BZ`` on Q2
//	flops +=  4 * m * (n-k) * k
//		- 2 * (n-k) * k * k
//		- 4 * m * (n-k)
//		- 3 * (n-k)
//		+ 5 * (n-k) * k ;
//
//      this is the LARF to create Q2
//	flops += 4 * (m-k+1)  * (n-k) - (n-k) ;
 
//      this is the cost for creating Q1
//	flops = ( 
//	        + 6 * (m-k) * k
//	        + 4 * k * k
//		- 3 * (m-k)
//		- 1 )
//		* k
//		/ (( long int )  3 );

//      this is the cost for creating Q2
	flops =  ( 4 * m - 2 * k + 1 ) * k * ( n - k ) ;


	return flops;

}
