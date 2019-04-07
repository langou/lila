#include "qr3.h"

unsigned long int flops_larft( int m, int k ){

	unsigned long int u_m, u_k, flops;

	u_m = ( unsigned long int ) m;
	u_k = ( unsigned long int ) k;

	flops = ((( unsigned long int ) 3) * u_m * u_k * u_k - u_k * u_k * u_k ) / (( unsigned long int ) 3) ;

	return flops;

}
