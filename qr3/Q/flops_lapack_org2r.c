#include "qr3.h"

long int flops_lapack_org2r( int m, int n, int k ){

	long int flops;

	flops = (( long int ) 0 );
	
	flops += 12 * m * n * k ;
	flops -=  6 * m * k * k ;
	flops -=  6 * n * k * k ;
	flops +=  4 * k * k * k ;
	flops -=  3 * m * k ;
	flops +=  3 * n * k ;
	flops -=  1 * k ;

	flops = flops / 3 ;

	return flops;

}
