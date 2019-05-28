#include "flops.h"

long int flops_qr2_dorgqr_check( int int_m, int int_n ){

	long int flops;
	long int n1, n2;
	long int m, n;

	m = ( long int ) int_m;

	n = ( long int ) int_n;

	flops = (( long int ) 0 );

	if( n <= 1 ){

		flops += m;

	} else {

		n1 = n/2;
		n2 = n-n1;

		flops += flops_qr2_dorgqr_check( m-n1, n2 );

//		flops += n1 * n1 * n2; // we do not do these flops because we are doing a BZ
		flops += 2 * n1 * n2 * (m-n1);
		flops += (n1-1) * n1 * n2 ; // we do these flops but we want to make sure we have the good formula
		flops += n1*n2;
		flops += 2 * n1 * n2 * (m-n1);
		flops += n1 * n1 * n2;

		flops += flops_qr2_dorgqr_check( m, n1 );

	}

	return flops;

}
