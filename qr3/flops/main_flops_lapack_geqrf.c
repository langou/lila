#include "flops.h"

int main(int argc, char ** argv) {

	int i, m, n, nb;
	
    	m  = 6;
    	n  = 5;
	nb = 1;

	for(i = 1; i < argc; i++){
		if( strcmp( *(argv + i), "-m") == 0) {
			m = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-n") == 0) {
			n = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-nb") == 0) {
			nb = atoi( *(argv + i + 1) );
			i++;
		}
	}

//	printf("%5d %5d %5d %15ld %15ld %3ld\n", m, n, nb, flops_lapack_geqrf_check( m, n, nb ), flops_lapack_geqr2_check( m, n ),  flops_lapack_geqrf_check( m, n, nb )-flops_lapack_geqr2_check( m, n )) ;

	printf("%5d %5d %5d %15ld %15ld %15ld\n", m, n, nb,
		flops_lapack_geqrf( m, n, nb ),
		flops_lapack_geqrf_check( m, n, nb ),
		flops_lapack_geqrf( m, n, nb )-flops_lapack_geqrf_check( m, n, nb )) ;

	printf("%5d %5d %5d %15ld %15ld %15ld\n", m, n, nb,
		flops_lapack_geqrf_bef( m, n, nb ),
		flops_lapack_geqrf_bef_check( m, n, nb ),
		flops_lapack_geqrf_bef( m, n, nb )-flops_lapack_geqrf_bef_check( m, n, nb )) ;

	return 0;

}
