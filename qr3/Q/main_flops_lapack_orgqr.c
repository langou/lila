#include "qr3.h"

int main(int argc, char ** argv) {

	int i, m, n, k, nb;
	
    	m  = 6;
    	n  = 5;
    	k  = 4;
	nb = 2;

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
		if( strcmp( *(argv + i), "-nb") == 0) {
			nb = atoi( *(argv + i + 1) );
			i++;
		}
	}

	printf("%5d %5d %5d %5d %15ld %15ld\n", m, n, k, nb,
//		flops_org2r( m, n, k ),
//		flops_lapack_orgqr( m, n, k, nb ),
		flops_lapack_orgqr_check( m, n, k, nb ),
		flops_lapack_orgqr_from_first_org2r( m, n, k, nb ) );
//		flops_lapack_orgqr_from_last_org2r( m, n, k, nb ) );

	return 0;

}
