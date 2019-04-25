#include "qr3.h"

long int flops_lapack_geqrf_check( int m, int n, int nb ){

	int ml, nl, ib, k, i, j;
	long int flops;

	flops = (( long int ) 0 );
	
	j = 0;
	i = 0;
	k = n; 

	ib = nb; if( k - i - nb < 0 ) ib = k - i;

	ml = m;
	nl = n;

// 	Within the while, ib is nb
//	ml = m - j*nb
//	nl = n - j*nb
	while( k - i > nb ){

//		GEQR2
//		flops += flops_lapack_larfg( ml );
		flops += flops_lapack_geqr2( ml, ib );

//		LARFT
		flops += flops_larft( ml, ib );

//		LARFB
//		flops += flops_lapack_larf( ml, nl-ib );
		flops += flops_lapack_larfb( ml, nl-ib, ib );

		ml -= ib;		
		nl -= ib;		
		i  += ib;		
		j++;

		ib = nb; if( k - i - nb < 0 ) ib = k - i;
	
	}	

//	GEQR2 cleanup

	flops += flops_lapack_geqr2( ml, nl );


	return flops;

}
