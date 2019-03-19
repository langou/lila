#include "lila.h"

int lila_dormqrf_w03( int *lila_param, int m, int n, int k, int i, int j, int mt, double *A, int lda, double *T, int ldt, double *work, int lwork ){

	int vrtq;
	vrtq = lila_param[4];

	// vrtq == 0  ==> VRTQ
	// vrtq == 1  ==>  VRT
	if(( vrtq == 0 ) || ( vrtq == 1 ) ){
		
		double *Aii, *Aij, *Tii;
		int ml, vb, ldwork, ii, jj, not_done, kk;

		vb = mt - (i % mt );

		if ( vb > k ) vb = k;

		Aii = A + i        + i*lda;
		Aij = A + i        + j*lda;
		Tii = T + (i % mt) + i*ldt;

		ml = m - i;
		kk =     1;

		ldwork   = mt;
		not_done =  1;

		while( not_done == 1 ){

			for( jj = 0; jj < n; jj++ ){
				for( ii = 0; ii < vb; ii++ ){
					work[ ii + jj * ldwork ] = Aij[ ii + jj * lda  ];
				}
			}
	
			cblas_dtrmm( CblasColMajor, CblasLeft, CblasLower, CblasTrans, CblasUnit, vb, n, (1.0e+00), Aii, lda, work, ldwork );
			cblas_dgemm( CblasColMajor, CblasTrans, CblasNoTrans, vb, n, ml-vb, (1.0e+00), Aii+vb, lda, Aij+vb, lda, (1.0e+00), work, ldwork );
			cblas_dtrmm( CblasColMajor, CblasLeft, CblasUpper, CblasTrans, CblasNonUnit, vb, n, (1.0e+00), Tii, ldt, work, ldwork );
			cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, ml-vb, n, vb, (-1.0e+00), Aii+vb, lda, work, ldwork, (1.0e+00), Aij+vb, lda );
			cblas_dtrmm( CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, vb, n, (+1.0e+00), Aii, lda, work, ldwork );


			for( jj = 0; jj < n; jj++ ){
				for( ii = 0; ii < vb; ii++ ){
					Aij[ ii + jj * lda  ] -= work[ ii + jj * ldwork ];
				}
			}

			if ( kk + vb - 1 == k ) {

				not_done = 0; 

			} else if ( kk + vb - 1 > k ) {

				printf("defensive programming, we should never have been there, abort\n"); return 0;

			} else {

				if( kk == 1 ) Tii = Tii - ( i%mt ) + vb * ldt;

				ml -= vb;
				kk += vb;

				Aii += vb * ( lda + 1 );
				Aij += vb;
				if(kk != 1+vb) Tii = Tii + vb * ldt;

				if ( ( kk + mt - 1 ) <= k ) vb = mt; else vb = k - kk + 1;

			}

		}	

	}

	return 0;

}
