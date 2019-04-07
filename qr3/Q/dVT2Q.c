#include "qr3.h"

int dVT2Q( int m, int n, double *Q, int ldq  ){

	int i, info;

//	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'L', m-1, n, A+1, lda, Q+1, ldq );
//	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'U', n, n, T, ldt, Q, ldq );

	info = ULTinU( n, Q, ldq, Q, ldq );

	cblas_dtrmm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, m-n, n, (-1.0e+00), Q, ldq, Q+n, ldq );

	info = mLUinA( n, Q, ldq );

	for(i = 0; i < n; i++) Q[ i + ldq * i ] = 1.00e+00 + Q[ i + ldq * i ];

	return 0;

}