#include "qr3.h"

long int flops_lapack_geqrf( int m, int n, int nb ){

	int ml, nl, ib, k, i, j;
	long int flops;

	flops = (( long int ) 0 );
	
	j = 0;
	i = 0;
	k = n; 

	ib = nb; if( k - i - nb < 0 ) ib = k - i;

	ml = m;
	nl = n;

	while( k - i > nb ){

//		GEQR2
//		flops += flops_lapack_geqr2( ml, ib );
//		flops += flops_lapack_geqr2_check( ml, ib );

//		LARFT
		flops += flops_larft( ml, ib );

//		LARFB
//		flops += flops_lapack_larfb( ml, nl-ib, ib );

//	flops += ib * ib * (nl-ib) ;            // 1: TRMM
//	flops += 2 * ib * (nl-ib) * ( ml-ib );   // 2: GEMM
	flops += (ib-1) * ib * (nl-ib) ;        // 3: TRMM (extra from a bunch of LARF)
//	flops += ib * (nl-ib) ;                // 3: TRMM
//	flops += 2 * ib * (nl-ib) * ( ml-ib );   // 4: GEMM
//	flops += ib * ib * (nl-ib);             // 5: TRMM


		ml -= ib;		
		nl -= ib;		
		i  += ib;		
		j++;

		ib = nb; if( k - i - nb < 0 ) ib = k - i;
	
	}	

//	GEQR2 cleanup

//	flops += flops_lapack_geqr2( ml, nl );
//	flops += flops_lapack_geqr2_check( ml, nl );

	flops += flops_lapack_geqr2( m, n);

	return flops;

}
