#include "qr3.h"

int main(int argc, char ** argv) {

	int i, m, n, k;
	
	m = 10;
    	n = 5;
	k = 5;

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

	printf("%d %d %d %lu\n", m, n, k, flops_org2r( m, n, k ) );

	return 0;

}
