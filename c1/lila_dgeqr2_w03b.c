#include "lila.h"

int lila_dgeqr2_w03b( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *TTT, int llldddttt, double *Q, int ldq, double *work, int lwork ){

	int info ; 
	double *tau=NULL;
	double *Aii, *Tii, *Qii;
	int ml;

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

/////////  This should create the good TTT now

	info = lila_dlarft_w03( m, n, i, mt, A, lda, TTT, llldddttt, tau);

//////////

//  	info = dgeqr3( ml, n, Aii, lda, Tii, ldt );
//	for(j = 0; j < n; j++) tau[j] = Tii[j+j*ldt];

	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml, n, Aii, lda, Qii, ldq );

	info = LAPACKE_dorgqr_work( LAPACK_COL_MAJOR, ml, n, n, Qii, ldq, tau, work, lwork );

// // // // //

//	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml, n, Asave, ml, Qii, ldq );

//	cblas_dsyrk( CblasColMajor, CblasUpper, CblasTrans, n, ml, 1.0e+00, Qii, ldq, 0e+00, R, n );

//	info = LAPACKE_dpotrf( LAPACK_COL_MAJOR, 'U', n, R, n ); 

//	cblas_dtrsm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, ml, n, 1.0e+00, R, n, Qii, lda );

/*
	printf("R second = \n");
	for( i1 = 0; i1 < n ; i1++){
	for( j1 = 0; j1 < n ; j1++){
		printf("%e,",R[ i1 + j1*ldt]);
	}
	printf("\n");
	}

	for( i1 = 0; i1 < n ; i1++){
		if ( R[ i1 + i1*n ] * Aii[ i1 + i1*lda ] < 0 ){
			for( j1 = 0; j1 < n ; j1++) R[ i1 + j1*n ] = - R[ i1 + j1*n ];
			for( j1 = 0; j1 < ml ; j1++) Qii[ j1 + i1*ldq ] = - Qii[ j1 + i1*ldq ];
		}
	}
*/


//	for( i1 = 0; i1 < n ; i1++){
//		if ( fabs(1 - Aii[ i1 + i1*lda ]) < fabs(1 + Aii[ i1 + i1*lda ]) ){
//			for( j1 = 0; j1 < n ; j1++) R[ i1 + j1*n ] = - R[ i1 + j1*n ];
//			for( j1 = 0; j1 < ml ; j1++) Qii[ j1 + i1*ldq ] = - Qii[ j1 + i1*ldq ];
//		}
		// This is the LU putting it right into A, L is lower unit and U is triu(A)
//		Aii[ i1 + i1*lda ] += 1; 
//		for( j1 = 0; j1 < m-i-i1; j1++ ) Aii[ j1 + i1*lda ] = Aii[ j1 + i1*lda ] / Aii[ i1 + i1*lda ]; 
//		cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, m-i-i1, n-i1, 1, (-1.0e+00), Aii+i1, lda, Aii+i1*lda, lda, (1.0e+00), Aii, lda ); 
//	}

////////////
//		This portion is what I added thinking this is how to apply 
//		T(itlo:ithi,jtlo:jthi) = triu( A(jtlo:jthi,jtlo:jthi) ) / ( (eye(vb,vb) + tril(A(jtlo:jthi,jtlo:jthi),-1) )');
//	
//	cblas_dtrsm( CblasColMajor, CblasRight, CblasLower, CblasTrans, CblasUnit, n, n, 1.0e+00, Aii, lda, Tii, ldt );
//
//	for( i1 = 0; i1 < n ; i1++){
//	for( j1 = i1; j1 < n ; j1++){
//		Aii[ i1 + j1*ldt ]= R[i1+j1*n];
//	}
//	}
//////////////
// this is what was there before I added the above section




//	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'U', n, n, R, n, Qii, lda );


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
