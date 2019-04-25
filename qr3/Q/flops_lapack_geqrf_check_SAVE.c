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
//		flops += flops_lapack_geqr2( ml, ib );
		//flops += 6 * ml * ib * ib - 2 * ib * ib * ib + 3 * ml * ib + 3 * ib * ib + 14 * ib;
		//flops += 6 * ml * nb * nb - 2 * nb * nb * nb + 3 * ml * nb + 3 * nb * nb + 14 * nb; 

		// ====>
		//flops += 
		//6 * ml * nb * nb 
		//- 2 * nb * nb * nb 
		//+ 3 * ml * nb 
		//+ 3 * nb * nb 
		//+ 14 * nb;

		// ====>
		flops += 
		6 * m * nb * nb 
		- 6 * j * nb * nb * nb
		- 2 * nb * nb * nb 
		+ 3 * m * nb 
		- 3 * j * nb * nb
		+ 3 * nb * nb 
		+ 14 * nb;

//		LARFT
//		flops += flops_lapack_larf( ml, ib );
		//flops += 4 * ml * ib + ib;
		//flops += 4 * ml * nb + nb;

		// ====>
		//flops += 
		//4 * ml * nb 
		//+ nb;

		// ====>
		flops += 
		4 * m * nb 
		- 4 * j * nb * nb
		+ nb;

//		LARFB
//		flops += flops_lapack_larfb( ml, nl-ib, ib );
		//flops += 4 * ml * ( nl - ib ) * ib - ( nl - ib ) * ib * ib + 2 * ( nl - ib ) * ib;
		//flops += 4 * ml * ( nl - nb ) * nb - ( nl - nb ) * nb * nb + 2 * ( nl - nb ) * nb;

		// ====>
		//flops += 
		//4 * ml * nl * nb 
		//- 4 * ml * nb * nb 
		//- nl * nb * nb 
		//+ nb * nb * nb 
		//+ 2 * nl * nb 
		//- 2 * nb * nb;

		// ====>
		//flops += 
		//4 * ( m - j*nb ) * ( n - j*nb ) * nb 
		//- 4 * ( m - j*nb ) * nb * nb 
		//- ( n - j*nb ) * nb * nb 
		//+ nb * nb * nb 
		//+ 2 * ( n - j*nb ) * nb
 
		// ====>
		flops += 
		4 * m * n * nb 
		- 4 * m * j*nb * nb 
		- 4 * j * nb * n * nb 
		+ 4 * j * nb * j * nb * nb 
		- 4 * m * nb * nb 
 		+ 4 * j * nb * nb * nb 
		- n * nb * nb 
		+ j * nb * nb * nb 
		+ nb * nb * nb 
		+ 2 * n * nb
		- 2 * j*nb * nb;

		ml -= ib;		
		nl -= ib;		
		i  += ib;		
		j++;

		ib = nb; if( k - i - nb < 0 ) ib = k - i;
	
	}	

//	GEQR2 cleanup
//	flops += flops_lapack_geqr2( ml, nl );
	//flops += 6 * ml * nl * nl - 2 * nl * nl * nl + 3 * ml * nl + 3 * nl * nl + 14 * nl;

	// ====>
	flops += 
	6 * ml * nl * nl 
	- 2 * nl * nl * nl 
	+ 3 * ml * nl 
	+ 3 * nl * nl 
	+ 14 * nl;

	printf(" j = %d, n %% nb = %d, \n", j, nb%n );

	return flops;

}
