#include "qr3.h"

unsigned long int flops_gemm( int m, int n, int k ){

	unsigned long int u_m, u_n, u_k, flops;

	u_m = ( unsigned long int ) m;
	u_n = ( unsigned long int ) n;
	u_k = ( unsigned long int ) k;

	flops = (( unsigned long int ) 2) * u_m * u_n * u_k;

	return flops;

}


