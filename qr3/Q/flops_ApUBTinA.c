#include "qr3.h"

unsigned long int flops_ApUBTinA( int m, int n ){

	unsigned long int u_m, u_n, flops;

	u_m = ( unsigned long int ) m;

	u_n = ( unsigned long int ) n;

	flops = u_m * u_m * u_n + u_m * u_n ;

	return flops;

}
