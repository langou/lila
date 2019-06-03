#include "../flops/flops.h"

int main(int argc, char ** argv) {

	int i, m, n, k;
	
    	m  = 6;
    	n  = 5;
    	k  = 4;

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

	printf("%5d %5d %5d %15ld %15ld ( %3ld )\n", m, n, k, flops_lapack_larfb( m, n, k ), flops_lapack_larfb_check( m, n, k ), flops_lapack_larfb( m, n, k ) - flops_lapack_larfb_check( m, n, k ) );

	return 0;

}
