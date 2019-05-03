#include "qr3.h"

long int flops_lapack_orgqr_from_first_org2r( int int_m, int int_n, int int_k, int int_b ){

	long int m, n, k, b, flops;
	long int kb;

	flops = (( long int ) 0 );

	m = ( long int ) int_m;
	n = ( long int ) int_n;
	k = ( long int ) int_k;
	b = ( long int ) int_b;

	if( k%b == 0 ) kb = k/b - 1; else kb = k/b;

//	flops += flops_org2r( m-k+ib, n-k+ib, ib );

	flops = ( (( long int ) 12 ) * (m-b*kb) * (n-b*kb) * (k-kb*b)
		- (( long int )  6 ) * (m-b*kb) * (k-kb*b) * (k-kb*b) 
		- (( long int )  6 ) * (n-b*kb) * (k-kb*b) * (k-kb*b) 
		+ (( long int )  4 ) * (k-kb*b) * (k-kb*b) * (k-kb*b)
		+ (( long int )  6 ) * (m-b*kb) * (k-kb*b) 
		- (( long int )  3 ) * (k-kb*b) * (k-kb*b) 
		- (( long int )  1 ) * (k-kb*b) 
		- (( long int )  3 ) ) 
		/ (( long int )  3 );

//	flops = ( (( long int ) 12 ) * (m-b*kb) * (n-b*kb) * (k-kb*b)
//		- (( long int )  6 ) * (m-b*kb) * (k-kb*b) * (k-kb*b) 
//		- (( long int )  6 ) * (n-b*kb) * (k-kb*b) * (k-kb*b) 
//	
//	
//		+ (( long int )   4 ) * k * k * k
//		- (( long int )  12 ) * k * k * kb * b
//		+ (( long int )  12 ) * k * kb * kb * b * b
//		- (( long int )   4 ) * kb * kb * kb * b * b * b
//	
//		+ (( long int )   6 ) * m * k 
//		- (( long int )   6 ) * m * kb * b 
//	
//		- (( long int )   6 ) * k * kb * b 
//		+ (( long int )   6 ) * kb * kb * b * b
//	
//		- (( long int )   3 ) * k * k 
//		+ (( long int )   6 ) * k * kb *b 
//		- (( long int )   3 ) * kb * kb * b * b 
//	
//		- (( long int )   1 ) * k 
//		+ (( long int )   1 ) * kb * b 
//		- (( long int )   3 ) ) 
//		/ (( long int )   3 );



	return flops;

}
