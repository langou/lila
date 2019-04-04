#include "qr3.h"

unsigned long int flops_trmm( int m, int n, char S ){

	unsigned long int u_m, u_n, flops;

	u_m = ( unsigned long int ) m;
	u_n = ( unsigned long int ) n;

	if( S == 'L' ){
		flops = (+2.0e00)*u_n*u_m*u_m;
	}

	if( S == 'R' ){
		flops = (+2.0e00)*u_n*u_n*u_m;
	}

	return flops;

}


