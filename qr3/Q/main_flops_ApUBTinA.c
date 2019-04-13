#include "qr3.h"

int main(int argc, char ** argv) {

	int i, m, n;
	
    	m = 6;
    	n = 5;

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

	printf("%d %d %lu %lu\n", m, n, flops_ApUBTinA( m, n ), flops_ApUBTinA_check( m, n ) );

	return 0;

}
