#include "../flops/flops.h"

int main(int argc, char ** argv) {

	int i, n;
	
    	n = 5;

	for(i = 1; i < argc; i++){
		if( strcmp( *(argv + i), "-n") == 0) {
			n = atoi( *(argv + i + 1) );
			i++;
		}
	}

	printf("%5d %15lu %15lu ( %3ld )\n", n, flops_V2N( n ), flops_V2N_check( n ), flops_V2N( n ) - flops_V2N_check( n ) );

	return 0;

}
