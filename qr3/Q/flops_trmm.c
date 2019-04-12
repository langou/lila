#include "qr3.h"

unsigned long int flops_trmm( char S, int m, int n ){

	unsigned long int u_m, u_n, flops;

	u_m = ( unsigned long int ) m;
	u_n = ( unsigned long int ) n;

	if( S == 'L' ){
		flops = u_n * u_m * u_m;
	}

	if( S == 'R' ){
		flops = u_n * u_n * u_m;
	}

	return flops;

}
