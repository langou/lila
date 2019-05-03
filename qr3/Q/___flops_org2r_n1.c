#include "qr3.h"

long int flops_org2r_n1( int int_m ){

	long int m, flops;

	m = ( long int ) int_m;

	flops = (( long int )  4 ) * m - (( long int )  3 ) ;

	return flops;

}
