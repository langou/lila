#include "qr3.h"

unsigned long int flops_N2T( int n ){

	unsigned long int u_n, flops;

	u_n = ( unsigned long int ) n;

	flops = ( u_n * u_n * u_n - u_n ) / (( unsigned long int ) 3 );

	return flops;

}
