#include "qr3.h"

long int flops_legacy_lapack_geqr2_check( int m, int n ){

	int ml, nl;
	long int flops;

	flops = (( long int ) 0 );
	
	ml = m;
	nl = n;

	while( 1 < nl  ){

//		GEQR2
		flops += flops_lapack_larfg( ml );
//		flops += flops_legacy_lapack_geqr2( ml, 1 );

//		LARF
		flops += flops_legacy_lapack_larf( ml, nl-1 );

		ml --;		
		nl --;		

	}	

//	GEQR2 cleanup

	flops += flops_lapack_larfg( ml );
//	flops += flops_legacy_lapack_geqr2( ml, 1 );

	return flops;

}
