#include "flops.h"

long int flops_dlarft3( int int_m, int int_n ){

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

	flops += flops_dlarft3( m, n1 );
	flops += flops_dlarft3( m-n1, n2 );

	flops += n1 * n2 * n2;

	flops += (( long int ) 2) * n1 * n2 * m;

	flops -= (( long int ) 2) * n1 * n2 * n;

	flops += n1 * n1 * n2;

	flops += n1 * n2 * n2;

	}	

	return flops;

}
