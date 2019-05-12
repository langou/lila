#include "qr3.h"

long int flops_lapack_org2r_Q1_check( int m, int n, int k ){

	int ml, nl;

	long int flops;

	if (n==1) return m;

	flops = (( long int ) 0 );
	
	ml = m-k+1;
	nl = n-k+1;

////////////////

//	flops += flops_lapack_larfb( ml, nl-1, 1 );
//	flops += flops_lapack_org2r_check( ml, 1, 1 );

//	flops += 4 * ml * (nl-1) - (nl-1) ;
	flops += ml ;

////////////////

	while( nl < n ){

		ml ++;		
		nl ++;
	
////////////////

//		flops += flops_lapack_larfb( ml, nl-1, 1 );

		flops += 4 * ml * (nl-1-(n-k)) - (nl-1-(n-k)) ;

////////////////

////////////////

//		flops += flops_lapack_org2r_check( ml, 1, 1 );

		flops += ml ;

////////////////

	}

	return flops;

}
