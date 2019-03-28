#include "lila.h"

int nV2T( int m, int n, double *A, int lda, double *tau, double *T, int ldt, double *work, int lwork ){

	double *Aii, *Tii, *taui;
	double *T0i, *Ai0, *Aj0, *Aji;
	int info, n1, n2, ii, jj; 

	if ( n <= 1 ) {

		T[0] = tau[0];

	} else {

		n1 = n / 2;
		n2 = n - n1;

		Aii = A + n1 + n1*lda;
		Tii = T + n1 + n1*ldt;
		taui = tau + n1;

		T0i = T + n1*ldt;
		Ai0 = A + n1;
		Aj0 = A + n1 + n2;
		Aji = A + n1 + n2 + n1*lda;

		info = nV2T( m, n1, A, lda, tau, T, ldt, work, lwork );
		info = nV2T( m-n1, n2, Aii, lda, taui, Tii, ldt, work, lwork );

		for( jj = 0; jj < n2; jj++ ){
			for( ii = 0; ii < n1; ii++ ){
				T0i[ ii + jj * ldt ] = Ai0[  jj + ii * lda  ];
			}
		}
		cblas_dtrmm( CblasColMajor, CblasRight, CblasLower, CblasNoTrans, CblasUnit, n1, n2, (+1.0e+00), Aii, lda, T0i, ldt );
		cblas_dgemm( CblasColMajor, CblasTrans, CblasNoTrans, n1, n2, m-n1-n2, (+1.0e+00), Aj0, lda, Aji, lda, (+1.0e+00), T0i, ldt );


		cblas_dtrmm( CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, n1, n2, (-1.0e+00), T, ldt, T0i, ldt );
		cblas_dtrmm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, n1, n2, (+1.0e+00), Tii, ldt, T0i, ldt );

	}

	return 0;

}


