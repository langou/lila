#include "qr2.h"

int qr2_dorgqr3_VS2Q( int m, int n, double *Q, int ldq, double *T, int ldt, double *tau ){


	double *Q11, *Q21;
	double *Q12, *Q22;
	double *T11, *T22;
	double *tau1, *tau2;

	int n1, n2;

	if( n <= 1 ){

//		LAPACKE_dorgqr_work( LAPACK_COL_MAJOR, m, 1, 1, Q, ldq, tau, work, lwork );

//		if ( m > 1 ) cblas_dscal( m-1, -(*tau), Q+1, 1);

		for( int i=1; i<m; i++)  Q[i] *= -(*tau);
		(*Q) = (+1.0e+00) - (*tau);


	} else {

		n1 = n/2;
		n2 = n-n1;

		Q11 = Q;
		Q21 = Q+n1;

		Q12 = Q+n1*ldq;
		Q22 = Q+n1+n1*ldq;

		T11 = T;
		T22 = T+n1+n1*ldt;

		tau1 = tau;
		tau2 = tau+n1;

		qr2_dorgqr3_VS2Q( m-n1, n2, Q22, ldq, T22, ldt, tau2 );

		cblas_dgemm( CblasColMajor, CblasTrans, CblasNoTrans, n1, n2, m-n1, (+1.0e+00), Q21, ldq, Q22, ldq, (+0.0e+00), Q12, ldq );
		cblas_dtrsm( CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, n1, n2, (+1.0e+00), T11, ldt, Q12, ldq );
		cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, m-n1, n2, n1, (-1.0e+00), Q21, ldq, Q12, ldq, (+1.0e+00), Q22, ldq );
		cblas_dtrmm( CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, n1, n2, (-1.0e+00), Q11, ldq, Q12, ldq );

		if (T11 != Q11 ) LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'U', n1, n1, T11, ldt, Q11, ldq );
		qr2_dVS2Q( m, n1, Q11, ldq );

	}

	return 0;

}
