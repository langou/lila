#include "lila.h"

unsigned long int flops_orgqr( int m, int n, int k ){

	unsigned long int u_m, u_n, u_k, flops;

	u_m = ( unsigned long int ) m;
	u_n = ( unsigned long int ) n;
	u_k = ( unsigned long int ) k;

	flops = (+4.0e00)*u_m*u_n*u_k - (+2.0e00)*(u_m+u_n)*u_k*u_k + (+4.0e00/3.0e00)*u_k*u_k*u_k + 2*u_m*u_k - u_k*u_k - (+1.0e00/3.0e00)*u_k;

	return flops;

}


