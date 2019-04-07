#include "qr3.h"

unsigned long int flops_syrk( int n, int k ){

	unsigned long int u_n, u_k, flops;

	u_n = ( unsigned long int ) n;
	u_k = ( unsigned long int ) k;

	flops = u_k * u_n * u_n;

	return flops;

}
