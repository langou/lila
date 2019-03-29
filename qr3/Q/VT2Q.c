#include "lila.h"

int VT2Q( int m, int n, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork ){

	int n1, n2, i, info; 
	double *Q11, *Q22, *Q02;

	if ( n <= 1 ) {

		info = LAPACKE_dlaset( LAPACK_COL_MAJOR, 'A', n, n, (0.0e+00), (0.0e+00), work, n);
		info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'U', n, n, Q, ldq, work, n );
		info = LAPACKE_dlaset( LAPACK_COL_MAJOR, 'U', n, n, (0.0e+00), (1.0e+00), Q, ldq);
		cblas_dtrmm( CblasColMajor, CblasRight, CblasLower, CblasTrans, CblasUnit, n, n, (+1.0e+00), Q, lda, work, n );
		cblas_dtrmm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, n, n, (-1.0e+00), work, n, Q, ldq );
	//	for(i = 0; i < n; i++) Q[ i + ldq * i ] = 1.00e+00 + Q[ i + ldq * i ];
		(*Q) = +1.00e+00 + (*Q);

	} else {

		n1 = n / 2;
		n2 = n - n1;

		Q11 = Q;
		Q02 = Q+n1*ldt;
		Q22 = Q+n1+n1*ldt;

		VT2Q( m, n1, A, lda, T, ldt, Q11, ldq, work, lwork );

		cblas_dtrmm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, m-n1, n1, (-1.0e+00), work, n1, Q11+n1, ldq );
	//	for(i = 0; i < n1; i++) Q11[ i + ldq * i ] = 1.00e+00 + Q11[ i + ldq * i ];

		VT2Q( m, n2, A, lda, T, ldt, Q22, ldq, work, lwork );

		cblas_dtrmm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, m-n2-n1, n2, (-1.0e+00), work, n2, Q22+n2, ldq );
	//	for(i = 0; i < n2; i++) Q22[ i + ldq * i ] = 1.00e+00 + Q22[ i + ldq * i ];

		// dormqrbz
		cblas_dgemm( CblasColMajor, CblasTrans, CblasNoTrans, n, n, m-n, (1.0e+00), A+n, lda, Q02+n, ldq, (0.0e+00), work, n );
		cblas_dtrmm( CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, n, n, (1.0e+00), T, ldt, work, n );
		LAPACKE_dlacpy( LAPACK_COL_MAJOR, 'A', n, n, work, n, Q02, ldq );
		cblas_dtrmm( CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, n, n, (-1.0e+00), A, lda, Q02, ldq );
		cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, m-n, n, n, (-1.0e+00), A+n, lda, work, n, (1.0e+00), Q02+n, ldq );



	}

	return 0;

}


