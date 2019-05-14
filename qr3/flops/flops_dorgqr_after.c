#include "flops.h"

long int flops_dorgqr_after( int m, int n, int k ){

	long int flops;

	flops = (( long int ) 0 );

	// copy V2^T into Q1
	//for ( i = 0; i < k; i++ ) for ( j = 0; j < n-k; j++ ) Q1[i+j*ldq] = A2[j+i*lda];
//	flops += n * k - k * k;

	// - T * V2^T -- note putting the minus sign on this term
	//cblas_dtrmm( CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit,  k, n-k, (-1.0e+00), T, ldt, Q1, ldq );
//	flops += n * k * k - k * k * k;
 
	// V2 * ( T * V2^T )
	//cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, n-k, n-k, k, (+1.0e+00), A2, lda, Q1, ldq, (+0.0e+00), Q2, ldq);
	//flops += flops_gemm( n-k, n-k, k );
//	flops += (n-k)*(n-k)*k;
//	flops += (n*n-n*k-k*n+k*k)*k;
//	flops += 2 * ( n * n * k - n * k * k - n * k * k + k * k * k );

	// Adding the identity block to Q2
	//for( j = 0; j < n-k; j++ ){ Q2[ j + j*ldq ] += (+1.0e+00); }
//	flops += n - k;

	// V3 * ( T * V2^T )
	//cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, m-n, n-k, k, (+1.0e+00), A3, lda, Q1, ldq, (+0.0e+00), Q3, ldq);
	//flops += flops_gemm( m-n, n-k, k );
//	flops += (m-n)*(n-k)*k;
//	flops += (m*n-m*k-n*n+n*k)*k;
//	flops += 2 * ( m * n * k - m * k * k - n * n * k + n * k * k );

	// V1 * ( T * V2^T ) 
	//cblas_dtrmm( CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, k, n-k, (+1.0e+00), A1, lda, Q1, ldq );
//	flops += n * k * k - k * k * k;


	//flops += 2 * m * n * k;
	//flops -= 2 * m * k * k;
	//flops += n * k;
	//flops -= k * k;
	//flops += n;
	//flops -= k;

	flops = (( long int ) 2 ) * ( m * n * k  - m * k * k ) + n * k - k * k + n - k;

	return flops;

}
