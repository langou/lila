#include "qr3.h"

unsigned long int flops_VT2Q( int m, int n ){

	unsigned long int u_m, u_n, flops;

	u_m = ( unsigned long int ) m;
	u_n = ( unsigned long int ) n;

	flops = u_m * u_n * u_n + ( u_n * u_n - u_n ) / (( unsigned long int ) 2 ) ;

	return flops;

}
