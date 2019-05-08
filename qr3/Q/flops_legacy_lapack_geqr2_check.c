#include "qr3.h"

long int flops_legacy_lapack_geqr2_check( int m, int n ){

	int ml, nl;
	long int flops;

	flops = (( long int ) 0 );
	
	ml = m;
	nl = n;

	while( 1 < nl  ){

////////////////

//		flops += flops_legacy_lapack_larfg( ml );

//		flops += flops_legacy_lapack_geqr2( ml, 1 );

		flops += 4 * ml + 5;

////////////////

////////////////

//		flops += flops_lapack_larfb( ml, nl-1, 1 );
//		flops += (nl-1);

//		flops += flops_legacy_lapack_larf( ml, nl-1 );
//		flops += (nl-1);

		flops += 4 * ml * (nl-1) ;

////////////////

		ml --;		
		nl --;		

	}	

////////////////

//	flops += flops_legacy_lapack_larfg( ml );

//	flops += flops_legacy_lapack_geqr2( ml, 1 );

	flops += 4 * ml + 5;

////////////////

	return flops;

}
