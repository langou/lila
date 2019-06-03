#include "flops.h"

long int flops_geqr3_noR( int int_m, int int_n ){

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

	flops += flops_geqr3_noR( m, n1 );

	//	for (i=0;i<n1;i++) for (j=0;j<n2;j++) T12[i+j*ldt] = A12[j*lda+i];
//		flops += n1 * n2;

	//	cblas_dtrmm ( CblasColMajor, CblasLeft, CblasLower, CblasTrans, CblasUnit, n1, n2, 1.0e+00, A11, lda, T12, ldt); 
		flops += flops_trmm( 'L', n1, n2 );

	//	cblas_dgemm( CblasColMajor, CblasTrans, CblasNoTrans, n1, n2, m-n1, +1.0e+00, A21, lda, A22, lda, +1.0e+00, T12, ldt);
		flops += flops_gemm( n1, n2, m-n1 );

	//	cblas_dtrmm ( CblasColMajor, CblasLeft, CblasUpper, CblasTrans, CblasNonUnit, n1, n2, 1.0e+00, T11, ldt, T12, ldt); 
		flops += flops_trmm( 'L', n1, n2 );

	//	cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, m-n1, n2, n1, -1.0e+00, A21, lda, T12, ldt, +1.0e+00, A22, lda);
		flops += flops_gemm( m-n1, n2, n1 );

	///////////////////////////////////////////////////
	//      These are what we save from not computing R
	//
	//	cblas_dtrmm ( CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, n1, n2, 1.0e+00, A11, lda, T12, ldt); 
	//	for (j=0;j<n2;j++) for (i=0;i<n1;i++) R12[j*ldr+i] = A12[j*lda+i] - T12[j*ldt+i]; 
	///////////////////////////////////////////////////


	flops += flops_geqr3_noR( m-n1, n2 );

	//	for (i=0;i<n1;i++) for (j=0;j<n2;j++) T12[i+j*ldt] = A21[j+i*lda];
//		flops += n1 * n2;

	//	cblas_dtrmm ( CblasColMajor, CblasRight, CblasLower, CblasNoTrans, CblasUnit, n1, n2, 1.0e+00, A22, lda, T12, ldt); 
		flops += flops_trmm( 'R', n1, n2 );

	//	cblas_dgemm( CblasColMajor, CblasTrans, CblasNoTrans, n1, n2, m-n, +1.0e+00, A21+n2, lda, A22+n2,lda, +1.0e+00, T12, ldt);
		flops += flops_gemm( n1, n2, m-n );

	//	cblas_dtrmm ( CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, n1, n2, -1.0e+00, T11, ldt, T12, ldt); 
		flops += flops_trmm( 'L', n1, n2 );

	//	cblas_dtrmm ( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, n1, n2, +1.0e+00, T22, ldt, T12, ldt); 
		flops += flops_trmm( 'R', n1, n2 );

	}


	return flops;

}
