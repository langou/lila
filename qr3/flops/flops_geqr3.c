#include "flops.h"

long int flops_geqr3( int m, int n ){

	long int flops;

	flops = (( long int ) 0 );
	
	// The extra flops to use the T matrix
	//flops = (  (n*n*n) - ((( long int ) 3 )*n*n) + ((( long int ) 2 )*n)  ) / (( long int ) 6 );
	// (1/6)*n*(n-1)*(n-2) ----- I'm finding this representation more stable. The above will result in negative values
	flops +=  ( n * ( n-(( long int ) 1 ) ) * ( n-(( long int ) 2 ) ) ) / (( long int ) 6 );


//	According to Elmroth and Gustavson:
//		flop_RGEQR3 = dlarfg + dlarft + ( the cost of doing the (k-1) update computations (Q^T * C) )


	flops += ( n * ( n-(( long int ) 1) ) * ( (( long int ) 6) * m // larft
	-( (( long int ) 2) * n 
	-(( long int ) 1) ) )  ) 
	/(( long int ) 6);

//	flops += (( long int ) 4) * m * n - n;   // larf
//	flops += (( long int ) 3) * m + 5;	 // larfg

	// lapack_dgeqr2
	flops += ( (( long int ) 6) * m * n * n 
	- (( long int ) 2) * n * n * n 
	+ (( long int ) 3) * m * n 
	+ (( long int ) 17) * n )
	/ (( long int ) 3);


	return flops;

}
