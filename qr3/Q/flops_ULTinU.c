#include "qr3.h"

unsigned long int flops_ULTinU( int n ){

	unsigned long int u_n, flops;

	u_n = ( unsigned long int ) n;

	flops = ( (( unsigned long int ) 2 ) * u_n * u_n * u_n + (( unsigned long int ) 3 ) * u_n * u_n - (( unsigned long int ) 5 ) * u_n ) / (( unsigned long int ) 6 );

	return flops;

}
