#include "lila.h"

int lila_dgeqr2_w02b( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork ){

	int info ; 
	double *tau=NULL;
	double *Aii, *Tii, *Qii;
	int ml;

	int i1, j1;

	tau = (double *) malloc( n * sizeof(double));

	Aii = A + i*lda + i;
	Qii = Q + i*ldq + i;
	Tii = T + i*ldt + i;

	ml = m - i;

	double *Asave, *R;

	Asave = (double *) malloc( ml * n * sizeof(double));
	R = (double *) malloc( n * n * sizeof(double));

//	int i1, j1;
//	for( i1 = 0; i1 < n ; i1++){
//	for( j1 = 0; j1 < n ; j1++){
//		R[ i1 + j1*n ]= -1.0e+00;
//	}
//	}

        cblas_dsyrk( CblasColMajor, CblasUpper, CblasTrans, n, ml, 1.0e+00, Aii, lda, 0e+00, R, n );

	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml, n, Aii, lda, Asave, ml );

  	info = LAPACKE_dgeqrf_work( LAPACK_COL_MAJOR, ml, n, Aii, lda, tau, work, lwork );

	info = LAPACKE_dlarft_work( LAPACK_COL_MAJOR, 'F', 'C', ml, n, Aii, lda, tau, Tii, ldt);

//  	info = dgeqr3( ml, n, Aii, lda, Tii, ldt );
//	for(j = 0; j < n; j++) tau[j] = Tii[j+j*ldt];

	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml, n, Aii, lda, Qii, ldq );

	info = LAPACKE_dorgqr_work( LAPACK_COL_MAJOR, ml, n, n, Qii, ldq, tau, work, lwork );

// // // // //

//	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml, n, Asave, ml, Qii, ldq );

//	cblas_dsyrk( CblasColMajor, CblasUpper, CblasTrans, n, ml, 1.0e+00, Qii, ldq, 0e+00, R, n );

//	info = LAPACKE_dpotrf( LAPACK_COL_MAJOR, 'U', n, R, n );

//	cblas_dtrsm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, ml, n, 1.0e+00, R, n, Qii, lda );

//	for( i1 = 0; i1 < n ; i1++){
//		if ( R[ i1 + i1*n ] * Aii[ i1 + i1*lda ] < 0 ){
//			for( j1 = 0; j1 < n ; j1++) R[ i1 + j1*n ] = - R[ i1 + j1*n ];
//			for( j1 = 0; j1 < ml ; j1++) Qii[ j1 + i1*ldq ] = - Qii[ j1 + i1*ldq ];
//		}
//	}

//	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'U', n, n, R, n, Aii, lda );


// continue here

// set T to zeros, copy R in T, trsm with Lower, Transpose, Unit, Right

//	for( i1 = 0; i1 < n ; i1++){
//	for( j1 = 0; j1 < n ; j1++){
//		Tii[ i1 + j1*ldt ]= 0.0e+00;
//	}
//	}

//	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'U', n, n, R, n, Tii, ldt );

//	cblas_dtrsm( CblasColMajor, CblasRight, CblasLower, CblasTrans, CblasUnit, n, n, 1.0e+00, Aii, lda, Tii, ldt );

	free( R );
	free( tau );
	free( Asave );

	return 0;

}
