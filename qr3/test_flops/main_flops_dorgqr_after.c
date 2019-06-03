#include "../flops/flops.h"

int main(int argc, char ** argv) {

	int i, m, n, k;
	
    	m  = 6;
    	n  = 5;
	k  = 3;

	for(i = 1; i < argc; i++){
		if( strcmp( *(argv + i), "-m") == 0) {
			m = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-n") == 0) {
			n = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-k") == 0) {
			k = atoi( *(argv + i + 1) );
			i++;
		}
	}

	printf("%5d %5d %5d %15ld %15ld (  %3ld )\n", m, n, k, flops_dorgqr_after( m, n, k ), flops_dorgqr_after_check( m, n, k ),  (  flops_dorgqr_after( m, n, k ) - flops_dorgqr_after_check( m, n, k ))) ;

	return 0;

}
