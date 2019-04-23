#include "qr3.h"

long int flops_org2r_check( int m, int n, int k ){

	long int flops;

	int k0, m1, n1, n2, ib;

	int nb = 1;

	flops = (( long int ) 0 );
	
	ib = k%nb; if (ib==0) ib=nb;

	k0 = k-ib;

	m1 = m-k0;
	n1 = n-k0;

	printf("\n");
//	flops += flops_org2r( m1, n1, ib );
	flops += 4 * m1 * (n1-1) + (n1-1);
	printf("larf(%3d,%3d)\n",m1,n1);
	flops += m1-1 ;
	flops ++ ;
	printf("org2r_n1(%d)\n",m1);

	while( k0 > 0 ){

		ib = nb; if( k0 - ib < 0 ) ib = k0;

		n2 = n1;
		m1 += ib;		
		n1 += ib;
		k0 -= ib;		

	//	flops += flops_larft( m1, ib );

		//flops += flops_lapack_larfb( m1, n2, ib );
		//printf(" %d %d %ld %ld %ld \n", m1 , n2, flops_lapack_larf( m1, n2 ), flops_lapack_larfb( m1, n2, ib ), flops_lapack_larf( m1, n2 ) - flops_lapack_larfb( m1, n2, ib ) );
		//flops += flops_lapack_larf( m1, n2 );

		flops += 4 * m1 * n2 + n2 ;
		printf("larf(%3d,%3d)\n",m1,n2);


		//flops += flops_org2r( m1, ib, ib );
//		printf(" %ld %ld \n ", flops_org2r_n1( m1 ), flops_org2r( m1, ib, ib ) );

		printf("org2r_n1(%d)\n",m1);
		flops += m1-1 ;
		flops ++ ;

// 		we would like that org2r returns larfg when nb=1
// 		we would like that larft returns 0 when nb=1
// 		we would like that larfb returns larf when nb=1

	}	

	return flops;

}
