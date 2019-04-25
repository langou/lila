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
//		flops +=   4 * nb * nb * nb * j * j; 


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
	flops += m * nb * nb * kb; 
	flops -= ( nb * nb * nb ) * ( ( kb * ( kb - (( long int ) 1 ) ) ) / (( long int ) 2 ) ); 
	flops -= m * nb * kb;
	flops += ( nb * nb ) * ( ( kb * ( kb - (( long int ) 1 ) ) ) / (( long int ) 2 ) ); 
	flops -= ( (( long int ) 2) / (( long int ) 6) )* nb * nb * nb * kb;
	flops += ( (( long int ) 1) / (( long int ) 6) ) * nb * nb * kb;
	flops -= ( ( (( long int ) 2) * nb * nb * nb - nb * nb ) / (( long int ) 6) ) * kb;
	flops += ( ( (( long int ) 2) * nb * nb - nb ) / (( long int ) 6) ) * kb;





//	LARFB Within the loop
	flops += (( long int ) 4 ) * m * n * nb * kb;
	flops -= (( long int ) 4 ) * m * nb * nb * ( ( kb * ( kb - (( long int ) 1 ) ) ) / (( long int ) 2 ) );
	flops -= (( long int ) 4 ) * n * nb * nb * ( ( kb * ( kb - (( long int ) 1 ) ) ) / (( long int ) 2 ) );
	flops += (( long int ) 4 ) * nb * nb * nb * ( ( kb * (kb - 1) * (2*kb - 1) ) / 6 );
	flops -= (( long int ) 4 ) * m * nb * nb * kb;
	flops += (( long int ) 4 ) * nb * nb * nb * ( ( kb * ( kb - (( long int ) 1 ) ) ) / (( long int ) 2 ) ) ;
	flops -= n * nb * nb * kb;
	flops += nb * nb * nb * ( ( kb * ( kb - (( long int ) 1 ) ) ) / (( long int ) 2 ) ) ;
	flops += nb * nb * nb * kb;
	flops += (( long int ) 2 ) * nb * n * kb;
	flops -= (( long int ) 2 ) * nb * nb * ( ( kb * ( kb - (( long int ) 1 ) ) ) / (( long int ) 2 ) );
	flops -= (( long int ) 2 ) * nb * nb * kb;




	return flops;

}
