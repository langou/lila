#include "qr3.h"

int main(int argc, char ** argv) {

	int i, n;
	
    	n = 5;

	for(i = 1; i < argc; i++){
		if( strcmp( *(argv + i), "-n") == 0) {
			n = atoi( *(argv + i + 1) );
			i++;
		}
	}

	printf("%d %lu\n", n, flops_V2N( n ) );

	return 0;

}
