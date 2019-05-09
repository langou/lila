#include "qr3.h"

long int flops_lapack_orgqr_bef( int m, int n, int k, int nb ){

	long int flops;

	int kb;

	if ( k%nb ) kb = k / nb; else kb = k / nb - 1; 

	flops = ( ( 6 * m + 6 * n - 6 * nb * kb - 2 * nb + 1 ) * ( nb * kb ) * ( nb - 1 ) ) / 6 ;

	return flops;

}
