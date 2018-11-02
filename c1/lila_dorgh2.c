#include "lila.h"

int lila_dorgh2( int m, int n, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork, int* S ){

//	this subroutine needs a workspace of size n-by-n
//	it is possible to reduce this workspace to nothing by changing the algorithm
//	but we think it is fine like this

	double *Qii, *Aii;
	int info, i1, j1, ml;

	//Aii = A + i*lda + i;
	//Qii = Q + i*ldq + i;
	//ml = m - i;

	Aii = A;
	Qii = Q;
	ml = m;

	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'U', n, n, Qii, ldq, T, ldt );

	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'U', n, n, T, ldt, work, n );
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'L', n-1, n, Aii+1, lda, work+1, n );

	for( i1 = 0; i1 < n ; i1++){

		if ( fabs(1 - work[ i1 + i1*n ]) < fabs(1 + work[ i1 + i1*n ]) ){

			S[ i1 ] = -1;

		} else {

			S[ i1 ] = 1;
			for( j1 = 0; j1 < n ; j1++) work[ j1 + i1*n ] = - work[ j1 + i1*n ];			
			for( j1 = n; j1 < ml ; j1++) Aii[ j1 + i1*lda ] = - Aii[ j1 + i1*lda ];			

		}

		work[ i1 + i1*n ] += 1.0e+00; 

		for( j1 = i1+1; j1 < n; j1++ ) work[ j1 + i1*n ] = work[ j1 + i1*n ] / work[ i1 + i1*n ];
		for( j1 = n; j1 < ml; j1++ ) Aii[ j1 + i1*lda ] = Aii[ j1 + i1*lda ] / work[ i1 + i1*n ];

		cblas_dger( CblasColMajor, n-i1-1, n-i1-1, (-1.0e+00), work + (i1+1) + i1*n, 1, work + i1 + (i1+1)*n, n, work + (i1+1) + (i1+1)*n, n);
		cblas_dger( CblasColMajor, ml-n, n-i1-1, (-1.0e+00), Aii + n + i1*lda, 1, work + i1 + (i1+1)*n, n, Aii + n + (i1+1)*lda, lda);

	}

	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'U', n, n, work, n, T, ldt ); 
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'L', n-1, n, work+1, n, Aii+1, lda );

	cblas_dtrsm( CblasColMajor, CblasRight, CblasLower, CblasTrans, CblasUnit, n, n, 1.0e+00, work, n, T, ldt );

	return 0;

}
