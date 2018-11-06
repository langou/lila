#include "lila.h"

int lila_dorgh2( int m, int n, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork, int* S ){

//	this subroutine needs a workspace of size n-by-n
//	it is possible to reduce this workspace to nothing by changing the algorithm
//	but we think it is fine like this

	int info, i1, j1;
int iiii;


//	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'U', n, n, Q, ldq, T, ldt );
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'U', n, n, T, ldt, work, n );
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'L', n-1, n, A+1, lda, work+1, n );

	for( i1 = 0; i1 < n ; i1++){

//for( j1 = 0; j1 < m; j1++){
//for( iiii = 0; iiii < n; iiii++){
//	printf(" %+9.4e ", A[ j1 + iiii*lda ]);
//}
//printf("\n");
//}
//printf("\n");



		//if ( fabs(1 - work[ i1 + i1*n ]) < fabs(1 + work[ i1 + i1*n ]) ){
//
			//S[ i1 ] = -1;
//
		//} else {

			S[ i1 ] = 1;
			for( j1 = 0; j1 < n ; j1++) work[ j1 + i1*n ] = - work[ j1 + i1*n ];			
			for( j1 = n; j1 < m ; j1++) A[ j1 + i1*lda ] = - A[ j1 + i1*lda ];			

		//}

		work[ i1 + i1*n ] += 1.0e+00; 

		for( j1 = i1+1; j1 < n; j1++ ) work[ j1 + i1*n ] = work[ j1 + i1*n ] / work[ i1 + i1*n ];
		for( j1 = n; j1 < m; j1++ ) A[ j1 + i1*lda ] = A[ j1 + i1*lda ] / work[ i1 + i1*n ];

		cblas_dger( CblasColMajor, n-i1-1, n-i1-1, (-1.0e+00), work + (i1+1) + i1*n, 1, work + i1 + (i1+1)*n, n, work + (i1+1) + (i1+1)*n, n);
		cblas_dger( CblasColMajor, m-n, n-i1-1, (-1.0e+00), A + n + i1*lda, 1, work + i1 + (i1+1)*n, n, A + n + (i1+1)*lda, lda);

	}

for( j1 = 0; j1 < n; j1++){
for( iiii = 0; iiii < n; iiii++){
	printf(" %+9.4e ", work[ j1 + iiii*n ]);
}
printf("\n");
}
printf("\n");



	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'U', n, n, work, n, T, ldt ); 
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'L', n-1, n, work+1, n, A+1, lda );

//for( j1 = 0; j1 < m; j1++){
//for( iiii = 0; iiii < n; iiii++){
//	printf(" %+9.4e ", A[ j1 + iiii*lda ]);
//}
//printf("\n");
//}
//printf("\n");

	cblas_dtrsm( CblasColMajor, CblasRight, CblasLower, CblasTrans, CblasUnit, n, n, 1.0e+00, work, n, T, ldt );

	return 0;

}
