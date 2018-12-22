#include "lila.h"

int lila_dgetrf_b( int m, int n, int i, int mt, double *A, int lda, double *Q, int ldq, double *work, int lwork ){

	double *Qii, *Aii, *U;
	int info, i1, j1, ml;

	ml = m - i;
	U = work;

	Aii = A + i*lda + i;
	Qii = Q + i*ldq + i;

	for( i1 = 0; i1 < n ; i1++){
		if ( fabs(1 - Aii[ i1 + i1*lda ]) < fabs(1 + Aii[ i1 + i1*lda ]) ){

			for( j1 = i1; j1 < n ; j1++) U[ i1 + j1*n ] = - U[ i1 + j1*n ];
			for( j1 = 0; j1 < ml ; j1++) Qii[ j1 + i1*ldq ] = - Qii[ j1 + i1*ldq ];

		} else {

			for( j1 = 0; j1 < ml ; j1++) Aii[ j1 + i1*lda ] = - Aii[ j1 + i1*lda ];			

		}

		Aii[ i1 + i1*lda ] += 1.0e+00; 
		for( j1 = i1+1; j1 < ml; j1++ ) Aii[ j1 + i1*lda ] = Aii[ j1 + i1*lda ] / Aii[ i1 + i1*lda ];
		cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, ml-i1-1, n-i1-1, 1, (-1.0e+00), Aii + (i1+1) + i1*lda, lda, Aii + i1 + (i1+1)*lda, lda, (1.0e+00), Aii + (i1+1) + (i1+1)*lda, lda ); 

	}

	return 0;

}
