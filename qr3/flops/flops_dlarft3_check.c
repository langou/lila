#include "flops.h"

long int flops_dlarft3_check( int int_m, int int_n ){

	long int flops;
	long int n1, n2;
	long int m, n;

	m = ( long int ) int_m;
	n = ( long int ) int_n;
	flops = (( long int ) 0 );

	if ( int_n <= 1){

	} else {

	n1 = int_n/2;
	n2 = int_n-n1;

	flops += flops_dlarft3_check( m, n1 );
	flops += flops_dlarft3_check( m-n1, n2 );

	//for (i=0;i<n1;i++) for (j=0;j<n2;j++) T12[i+j*ldt] = A21[j+i*lda];
//	flops += n1 * n2;

	//flops += flops_trmm( 'R', n1, n2 );
	flops += n1 * n2 * n2;
	//flops += flops_gemm( n1, n2, m-n );
	flops += (( long int ) 2) * n1 * n2 * (m - n);
	//flops += flops_trmm( 'L', n1, n2 );
	flops += n1 * n1 * n2;
	//flops += flops_trmm( 'R', n1, n2 );
	flops += n1 * n2 * n2;

	}	

	return flops;

}
