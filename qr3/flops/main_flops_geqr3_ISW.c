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

	printf("%5d %5d %15ld %15ld ( %15ld)\n", m, n, flops_geqr3_check(m,n), flops_geqr3(m,n), flops_geqr3_check(m,n) - flops_geqr3(m,n) ); // nothing to do with ISW - debugging

	printf("%5d %5d %15ld %15ld ( %15ld)\n", m, n, flops_geqr3_ISW_check(m,n), flops_geqr3_ISW(m,n), flops_geqr3_ISW_check(m,n) - flops_geqr3_ISW(m,n) ) ;

	printf("%5d %5d %15ld %15ld ( %15ld)\n", m, n, flops_geqr3_check(m,n), ( flops_geqr3_ISW(m,n) + flops_geqr3_ISW_save(m,n) ), flops_geqr3_check(m,n) - ( flops_geqr3_ISW(m,n) + flops_geqr3_ISW_save(m,n) ) ) ;

	printf("\n    qr3: %15ld   ISW: %15ld   ISW_flops_saved: %15ld\n", flops_geqr3_check(m,n), flops_geqr3_ISW_check(m,n), flops_geqr3_ISW_save(m,n) );
	printf("qr3-ISW: %15ld\n", flops_geqr3_check(m,n)-flops_geqr3_ISW_check(m,n) );

	return 0;

}
