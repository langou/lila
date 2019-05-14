#include "flops.h"

long int flops_lapack_geqrf_bef( int m, int n, int nb ){

	long int flops;

	int kb;

	flops = (( long int ) 0 );
	
	if ( n%nb ) kb = n / nb; else kb = n / nb - 1; 

	flops += m * ( nb * kb )  * ( nb - 1 );

	flops += ( n - nb * kb ) * ( nb * kb )  * ( nb - 1 ) ;

	flops -= (  ( 2 * nb - 1 ) * ( nb * kb ) * ( nb - 1 ) ) / 6;


	return flops;

}
