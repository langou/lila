#include "qr2.h"

int qr2_aux_dVT2Q( int m, int n, double *Q, int ldq  ){

	int i, info;

	info = qr2_aux_dULTinU( n, Q, ldq, Q, ldq );

	cblas_dtrmm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, m-n, n, (-1.0e+00), Q, ldq, Q+n, ldq );

	info = qr2_aux_dmLUinA( n, Q, ldq );

	for(i = 0; i < n; i++) Q[ i + ldq * i ] = 1.00e+00 + Q[ i + ldq * i ];

	return 0;

}
