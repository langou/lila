#include "../flops/flops.h"

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
		flops_qr2_dorgqr_check( m, n ),
		flops_lapack_org2r_check( m, n, n )-(n*n-n)/2,
		flops_qr2_dorgqr_check( m, n )-flops_lapack_org2r_check( m, n, n )+(n*n-n)/2);

//	the cost of org3r is 
//		org2r 
//		+ useT
//		- saveBZ
//	but since [ useT + n(n-1)/2  = saveBZ ]
//	we have that
//	the cost of org3r is 
//		org2r - n(n-1)/2


	printf("%5d %5d %15ld %15ld %15ld\n", m, n,
		flops_qr2_dorgqr_check( m, n ),
		flops_qr2_dorgqr( m, n ),
		flops_qr2_dorgqr_check( m, n )-flops_qr2_dorgqr( m, n ));


	return 0;

}
