#include "qr3.h"

int dV2Q( int m, int n, int k, int nb, int nbmin, int nx, double *A, int lda, double *tau, double *work, int lwork, int info ){

	double *Akk, *A0k, *Aik, *tauk;
	int kk, ml, nl, ib, i, j, flops;
	
	flops = 0;

	kk    = n;

	ml  = m-kk;
	nl  = n-kk;

	A0k  = A+kk*lda;
	Akk  = A+kk+kk*lda;
	tauk = tau+kk;

	ib = nb; if( kk - ib < 0 ) ib = n - nl;

	A0k  -= ib*lda;	
	Akk  -= ib*(1+lda);
	tauk -= ib;

	ml += ib;		
	nl += ib;		
	kk -= ib;

//	use unblocked code for the last or only block

//	flops += flops_org2r( ml, ib, ib );
//	dorg2r_( &ml, &ib, &ib, Akk, &lda, tauk, work, &lwork );

	for(i=0;i<ib;i++){ for(j=0;j<i;j++){ Akk[j+i*lda] = Akk[i+j*lda];}}

	flops += flops_V2N( ib );
	dV2N( ib, Akk, lda );

	flops += flops_syrk( ib, ml-ib );
	cblas_dsyrk( CblasColMajor, CblasUpper, CblasTrans, ib, ml-ib, (+1.0e+00), Akk+ib, lda, (+1.0e+00), Akk, lda );

	dN2T( ib, tauk, Akk, lda );

	dVT2Q( ml, ib, Akk, lda );

//	Use blocked code

	while( kk > 0 ){

		ib = nb; if( kk - ib < 0 ) ib = n - nl;

		A0k  -= ib*lda;
		Akk  -= ib*(1+lda);
		Aik   = Akk+ib*lda;
		tauk -= ib;

		ml += ib;		
		nl += ib;		
		kk -= ib;
	
//		Constructing T (in the upper part of the panel)

		for(i=0;i<ib;i++){ for(j=0;j<i;j++){ Akk[j+i*lda] = Akk[i+j*lda];}}

		flops += flops_V2N( ib );
		dV2N( ib, Akk, lda );

		flops += flops_syrk( ib, ml-ib );
		cblas_dsyrk( CblasColMajor, CblasUpper, CblasTrans, ib, ml-ib, (+1.0e+00), Akk+ib, lda, (+1.0e+00), Akk, lda );

		dN2T( ib, tauk, Akk, lda );

//	 	Popping Q up to current world - (bz)

//		Q1 = V2^T * Q2
		flops += flops_gemm( ib, nl-ib, nl-ib );
		cblas_dgemm( CblasColMajor, CblasTrans, CblasNoTrans, ib, nl-ib, nl-ib, (+1.0e+00), Akk+ib, lda, Aik+ib, lda, (+0.0e+00), Aik, lda );
//		Q1 = Q1 + V3^T * Q3
		flops += flops_gemm( ib, nl-ib, ml-nl );
		cblas_dgemm( CblasColMajor, CblasTrans, CblasNoTrans, ib, nl-ib, ml-nl, (+1.0e+00), Akk+nl, lda, Aik+nl, lda, (+1.0e+00), Aik, lda );
//		Q1 = T * Q1
		flops += flops_gemm( ib, nl-ib, 'L' );
		cblas_dtrmm( CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, ib, nl-ib, (+1.0e+00), Akk, lda, Aik, lda );
//		Q2 = Q2 - V2 * Q1
//		Q3   Q3   V3
		flops += flops_gemm( ml-ib, nl-ib, ib );
		cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, ml-ib, nl-ib, ib, (-1.0e+00), Akk+ib, lda, Aik, lda, (+1.0e+00), Aik+ib, lda );
//		Q1 = - V1 * Q1
		flops += flops_gemm( ib, nl-ib, 'L' );
		cblas_dtrmm( CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, ib, nl-ib, (-1.0e+00), Akk, lda, Aik, lda );

//		Constructing (ml x ib) panel of Q

//		flops += flops_org2r( ml, ib, ib );
//		dorg2r_( &ml, &ib, &ib, Akk, &lda, tauk, work, &lwork );
		dVT2Q( ml, ib, Akk, lda );
//
	}	

	printf("flops in total = %d\n", flops);
	return flops;

}
