#include "qr3.h"

unsigned long int flops_lapack_larfb( int m, int n, int k ){

	unsigned long int u_m, u_n, u_k, flops;

	u_m = ( unsigned long int ) m;
	u_n = ( unsigned long int ) n;
	u_k = ( unsigned long int ) k;

	flops = (( unsigned long int ) 4) * u_k * u_m * u_n - u_k * u_k * u_n + u_k * u_n ;

	return flops;

}
