#include "qr3.h"

long int flops_lapack_geqr2_check( int m, int n ){

	int ml, nl;
	long int flops;

	flops = (( long int ) 0 );
	
	ml = m;
	nl = n;

	while( 1 < nl  ){

//		GEQR2
//		flops += flops_lapack_larfg( ml );
		flops += flops_lapack_geqr2( ml, 1 );

//		LARFT
		flops += flops_larft( ml, 1 );

//		LARFB
//		flops += flops_lapack_larf( ml, nl-ib );
		flops += flops_lapack_larfb( ml, nl-1, 1 );

		ml --;		
		nl --;		

	}	

//	GEQR2 cleanup

//	flops += flops_lapack_larfg( ml );
	flops += flops_lapack_geqr2( ml, 1 );

	return flops;

}
