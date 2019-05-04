#include "qr3.h"

long int flops_legacy_lapack_org2r_check( int m, int n, int k ){

	long int flops;

	int m1, n1;
	int i;

	flops = (( long int ) 0 );
	
	m1 = m-k+1;
	n1 = n-k+1;

////////////////

//	flops += flops_legacy_lapack_org2r( m1, n1, 1 );

	flops += flops_legacy_lapack_larf( m1, n1-1 );
	flops += m1-1 ;
	flops ++ ;

//	flops += 4 * m1 * (n1-1) + (n1-1);
//	flops += m1-1 ;
//	flops ++ ;

////////////////

	for( i = 0; i < k-1; i++ ){

		m1 ++;		
		n1 ++;

////////////////

		flops += flops_legacy_lapack_larf( m1, (n1-1) );

//		flops += 4 * m1 * (n1-1) + (n1-1) ;

////////////////

////////////////

//		flops += flops_legacy_lapack_org2r( m1, 1, 1 );

		flops += m1-1 ;
		flops ++ ;

////////////////

	}	

	return flops;

}
