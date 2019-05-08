#include "qr3.h"

long int flops_legacy_lapack_org2r_check( int m, int n, int k ){

	int ml, nl;
	long int flops;

	flops = (( long int ) 0 );
	
	ml = m-k+1;
	nl = n-k+1;

////////////////

//	flops += flops_legacy_lapack_org2r( ml, nl, 1 );

//	flops += flops_legacy_lapack_larf( ml, nl-1 );
//	flops += ml-1 ;
//	flops ++ ;

	flops += 4 * ml * (nl-1) + (nl-1);
	flops += ml ;

////////////////

	while( ml < m  ){


		ml ++;		
		nl ++;

////////////////

//		flops += flops_lapack_larfb( ml, nl-1, 1 );
//		flops += 2*(nl-1) ;

//		flops += flops_legacy_lapack_larf( ml, (nl-1) );
//		flops += 2*(nl-1) ;

		flops += 4 * ml * (nl-1) + (nl-1) ;

////////////////

////////////////

//		flops += flops_legacy_lapack_org2r( ml, 1, 1 );

		flops += ml ;

////////////////

	}	

	return flops;

}
