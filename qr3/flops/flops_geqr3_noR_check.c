#include "flops.h"

long int flops_geqr3_noR_check( int int_m, int int_n ){

	long int flops;
	long int n1, n2;
	long int m, n;

	m = ( long int ) int_m;
	n = ( long int ) int_n;
	flops = (( long int ) 0 );

	if ( n == 1){

		flops += flops_lapack_larfg( m );

	} else {

	n1 = n/2;
	n2 = n-n1;

	flops += flops_geqr3_noR_check( m, n1 );

		flops += flops_trmm( 'L', n1, n2 );
		flops += flops_gemm( n1, n2, m-n1 );
		flops += flops_trmm( 'L', n1, n2 );
		flops += flops_gemm( m-n1, n2, n1 );

	///////////////////////////////////////////////////
	//      These are what we save from not computing R
	//
	//	cblas_dtrmm ( CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, n1, n2, 1.0e+00, A11, lda, T12, ldt); 
	//	for (j=0;j<n2;j++) for (i=0;i<n1;i++) R12[j*ldr+i] = A12[j*lda+i] - T12[j*ldt+i];
	// 
	///////////////////////////////////////////////////


	flops += flops_geqr3_noR_check( m-n1, n2 );

		flops += flops_trmm( 'R', n1, n2 );
		flops += flops_gemm( n1, n2, m-n );
		flops += flops_trmm( 'L', n1, n2 );
		flops += flops_trmm( 'R', n1, n2 );

	}


	return flops;

}
