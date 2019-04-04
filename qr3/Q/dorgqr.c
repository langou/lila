#include "qr3.h"

int dorgqr( int m, int n, int k, int nb, int nbmin, int nx, double *A, int lda, double *tau, double *work, int lwork, int info ){

	double *Akk, *A0k, *Aik, *T, *tauk;
	int kk, ldt, ml, nl, ib, i, j, flops;
	
	T     = work;
	ldt   = nb;
	work  = work+nb*(1+lwork);
	lwork = lwork-nb;

	flops = 0;
	kk    = n;

	ml  = m-kk;
	nl  = n-kk; // 0

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
//
//	Set A(1:kk,kk+1:n) to zero - this is part of LAPACK. Comment out to compare with our method
//
	//for( i = 0; i < kk; i++){ for( j = 0; j < ib; j++ ){ A0k[i+j*lda] = (+0.0e00); } }
//
//	use unblocked code for the last or only block
//
	flops += flops_orgqr( ml, ib, ib );
	info   = LAPACKE_dorgqr_work( LAPACK_COL_MAJOR, ml, ib, ib, Akk, lda, tauk, work, lwork );
//
//	Use blocked code
//
	while( kk > 0 ){

		ib = nb; if( kk - ib < 0 ) ib = n - nl;

		A0k  -= ib*lda;
		Akk  -= ib*(1+lda);
		Aik   = Akk+ib*lda;
		tauk -= ib;

		ml += ib;		
		nl += ib;		
		kk -= ib;
	
//	
//		LAPACK method for constructing T
//
	//	info = LAPACKE_dlarft_work( LAPACK_COL_MAJOR, 'F', 'C', ml, ib, Akk, lda, tauk, T, ldt);
//
//		Our 'faster' method for constructing T
//
		for(i=0;i<ib;i++){T[i+i*ldt] = 1.0e+00; for(j=0;j<i;j++){ T[j+i*ldt] = Akk[i+j*lda];}}
		xV2N( ib, Akk, lda );
		flops += flops_syrk( ib, ml-ib );
		cblas_dsyrk( CblasColMajor, CblasUpper, CblasTrans, ib, ml-ib, (+1.0e+00), Akk+ib, lda, (+1.0e+00), T, ldt );
		xN2T( ib, tauk, T, ldt );

		//info = LAPACKE_dlaset( LAPACK_COL_MAJOR, 'A', ib, nl-ib, (0.0e+00), (0.0e+00), Akk+ib*lda, lda);
	//	info = LAPACKE_dlarfb_work( LAPACK_COL_MAJOR, 'L', 'N', 'F', 'C', ml, nl-ib, ib, Akk, lda, T, ldt, Akk+ib*lda, lda, work, lwork);

//	 	Popping Q up to current world - (bz)
//
//		Q1 = V2^T * Q2
		flops += flops_gemm( ib, nl-ib, nl-ib );
		cblas_dgemm( CblasColMajor, CblasTrans, CblasNoTrans, ib, nl-ib, nl-ib, (+1.0e+00), Akk+ib, lda, Aik+ib, lda, (+0.0e+00), Aik, lda );
//		Q1 = Q1 + V3^T * Q3
		flops += flops_gemm( ib, nl-ib, ml-nl );
		cblas_dgemm( CblasColMajor, CblasTrans, CblasNoTrans, ib, nl-ib, ml-nl, (+1.0e+00), Akk+nl, lda, Aik+nl, lda, (+1.0e+00), Aik, lda );
//		Q1 = T * Q1
		flops += flops_gemm( ib, nl-ib, 'L' );
		cblas_dtrmm( CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, ib, nl-ib, (+1.0e+00), T, ldt, Aik, lda );
//		Q2 = Q2 - V2 * Q1
//		Q3   Q3   V3
		flops += flops_gemm( ml-ib, nl-ib, ib );
		cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, ml-ib, nl-ib, ib, (-1.0e+00), Akk+ib, lda, Aik, lda, (+1.0e+00), Aik+ib, lda );
//		Q1 = - V1 * Q1
		flops += flops_gemm( ib, nl-ib, 'L' );
		cblas_dtrmm( CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, ib, nl-ib, (-1.0e+00), Akk, lda, Aik, lda );
//
//		Constructing (ml x ib) panel of Q
//
		flops += flops_orgqr( ml, ib, ib );
		info   = LAPACKE_dorgqr_work( LAPACK_COL_MAJOR, ml, ib, ib, Akk, lda, tauk, work, lwork );
//
//		Set rows 1:kk of current block to zero - this is part of LAPACK. Comment out to compare with our method
//
	//	for( i = 0; i < kk; i++){ for( j = 0; j < ib; j++ ){ A0k[i+j*lda] = (+0.0e00); } }

	}	

	printf("flops in total = %d\n", flops);
	return flops;

}
