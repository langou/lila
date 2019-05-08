#include "qr3.h"

int main(int argc, char ** argv) {

	int i, m, n;
	
    	m  = 10;
    	n  =  6;

	for(i = 1; i < argc; i++){
		if( strcmp( *(argv + i), "-m") == 0) {
			m = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-n") == 0) {
			n = atoi( *(argv + i + 1) );
			i++;
		}
	}

	printf("%5d %5d %15ld %15ld %3ld\n", m, n,
			flops_legacy_lapack_geqr2( m, n ),
			flops_legacy_lapack_geqr2_check( m, n ), 
			flops_legacy_lapack_geqr2( m, n )-flops_legacy_lapack_geqr2_check( m, n )) ;

	return 0;

}
