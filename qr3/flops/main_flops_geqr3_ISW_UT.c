#include "flops.h"

int main(int argc, char ** argv) {

	int i, m, n;
	
    	m  = 6;
    	n  = 5;

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

	printf("%5d %5d %15ld %15ld %15ld\n", m, n,
		flops_geqr3_ISW_UT_check( m, n ),
		flops_geqr3_ISW_UT( m, n ),
		flops_geqr3_ISW_UT_check( m, n )- flops_geqr3_ISW_UT( m, n ));

	return 0;

}
