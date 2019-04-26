#include "qr3.h"

long int flops_lapack_geqrf( int m, int n, int nb ){

	int ml, nl, ib, k, i, j, kb;
	long int flops;

	flops = (( long int ) 0 );
	
	j = 0;
	i = 0;
	k = n; 

	ib = nb; if( k - i - nb < 0 ) ib = k - i;

	ml = m;
	nl = n;

	if( k%nb == 0 ) kb = k/nb - 1; else kb = k/nb;


// 	Within the while, ib is nb
//	ml = m - j*nb
//	nl = n - j*nb
	while( k - i > nb ){

//		GEQR2
//		flops += flops_lapack_larfg( ml );
//		flops += flops_lapack_geqr2( ml, ib );

//		LARFT
//		flops += flops_larft( ml, ib );


//		LARFB
//		flops += flops_lapack_larf( ml, nl-ib );
//		flops += flops_lapack_larfb( ml, nl-ib, ib );

		ml -= ib;		
		nl -= ib;		
		i  += ib;		
		j++;

		ib = nb; if( k - i - nb < 0 ) ib = k - i;
	
	}	


//	GEQR2 cleanup

//	flops += flops_lapack_geqr2( ml, nl );
	flops += ( (( long int ) 6) * ( m - j*nb ) * ( n - j*nb ) * ( n - j*nb ) 
	- (( long int ) 2) * ( n - j*nb ) * ( n - j*nb ) * ( n - j*nb ) 
	+ (( long int ) 3) * ( m - j*nb ) * ( n - j*nb ) 
	+ (( long int ) 3) * ( n - j*nb ) * ( n - j*nb ) 
	+ (( long int ) 14) * ( n - j*nb ) )
	/ (( long int ) 3);

//	GEQR2 within the while
	flops -= ( ( (( long int ) 2 ) * nb * nb * nb * kb ) / (( long int ) 3 ) );
	flops += nb * nb * kb;
	flops += ( ( (( long int ) 14 ) * nb * kb ) / (( long int ) 3 ) );
	flops += ( (( long int ) 3 ) * m * nb * kb ) / (( long int ) 3 );
	flops += ( (( long int ) 6 ) * m * nb * nb * kb ) / (( long int ) 3 );
	flops -= ( ( kb * ( kb - (( long int ) 1 ) ) ) / (( long int ) 2 ) ) * ( ( (( long int ) 6 ) * nb * nb * nb ) / (( long int ) 3 ) );  
	flops -= ( ( kb * ( kb - (( long int ) 1 ) ) ) / (( long int ) 2 ) ) * ( ( (( long int ) 3 ) * nb * nb ) / (( long int ) 3 ) );  

//	LARFT within the loop
	int flops_from_larft = 0;

	flops_from_larft += m * nb * nb * kb; 
	flops_from_larft -= ( nb * nb * nb ) * ( ( kb * ( kb - (( long int ) 1 ) ) ) / (( long int ) 2 ) ); 
	flops_from_larft -= m * nb * kb;
	flops_from_larft += ( nb * nb ) * ( ( kb * ( kb - (( long int ) 1 ) ) ) / (( long int ) 2 ) ); 
	flops_from_larft -= ( (( long int ) 2) / (( long int ) 6) )* nb * nb * nb * kb;
	flops_from_larft += ( (( long int ) 1) / (( long int ) 6) ) * nb * nb * kb;
	flops_from_larft -= ( ( (( long int ) 2) * nb * nb * nb - nb * nb ) / (( long int ) 6) ) * kb;
	flops_from_larft += ( ( (( long int ) 2) * nb * nb - nb ) / (( long int ) 6) ) * kb;

//	flops +=  flops_from_larft;


//	LARFB Within the loop

//	this is the flops in larfb due to the ( 3: TRMM )
	int flops_from_larfb_1 = 0;
//	flops_from_larfb_1 = ( 2 * n * nb * nb * kb - nb * nb * nb * kb * kb - nb * nb * nb * kb ) /2 ;
	flops_from_larfb_1 = ( n * n - n ) /2 ; // this the formula for when nb = 1
//	flops_from_larfb_1 += n * (n - kb * nb ) ; // this the formula for when nb = 1

	flops +=  flops_from_larfb_1;

//	this is the flops in larfb due to the (1,2,4,5,6) - everything but ( 3: TRMM )
	int flops_from_larfb_0 = 0;
	flops_from_larfb_0 += 12 * m * n * nb * kb ;               //   4 * m * n^2 
	flops_from_larfb_0 -=  6 * m * nb * nb * kb * kb;          // - 2 * m * n^2
	flops_from_larfb_0 -=  6 * m * nb * nb * kb;               // - 2 * m * n * b
	flops_from_larfb_0 +=  6 * n * nb * kb;                    // + 2 * n^2
	flops_from_larfb_0 -=  6 * n * nb * nb * kb * kb ;         // - 2 * n^3
	flops_from_larfb_0 +=  4 * nb * nb * nb * kb * kb * kb;    // + 4/3 * n^3
	flops_from_larfb_0 +=  3 * nb * nb * nb * kb * kb ;        // + 1 * n^2 * b
	flops_from_larfb_0 -=  1 * nb * nb * nb * kb ;             // - 1 * n * b^2
	flops_from_larfb_0 -=  3 * nb * nb * kb * kb;              // - 3 * n^2
	flops_from_larfb_0 -=  3 * nb * nb * kb;                   // - 1 * n * b

	flops +=  flops_from_larfb_0/3;

	return flops;

}
