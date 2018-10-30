#include "lila.h"

int lila_dgeqr2_w03b( int m, int n, int i, int mt, double *A, int lda, double *TTT, int llldddttt, double *Q, int ldq, double *work, int lwork ){

	int info ; 
	double *tau=NULL;
	double *Aii, *Qii, *TTTii, *R;
	int ml, i1, j1;

	R = (double *) malloc( n * n * sizeof(double));
	tau = (double *) malloc( n * sizeof(double));

	Aii = A + i*lda + i;
	Qii = Q + i*ldq + i;
	TTTii = TTT + i*llldddttt + i;
	
	ml = m - i;



//	double *Asave, *Qsave, *Asaveii;
	double *Asave;
//	int ldasave;
	Asave = (double *) malloc( ml * n * sizeof(double));
//	Qsave = (double *) malloc( ml * n * sizeof(double));
// // // // //
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml, n, Aii, lda, Asave, ml ); 
// // // // //


// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
// // // // //    This block is to make sure we still have a working code. I copy everything else to work with and will place pieces back to see if I can get it working.
// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //

//	we do this line to get tau
	info = LAPACKE_dgeqrf_work( LAPACK_COL_MAJOR, ml, n, Asave, ml, tau, work, lwork );
//	info = lila_dlarft_w03( m, n, i, mt, A, lda, TTT, llldddttt, tau);
//	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml, n, Aii, lda, Qii, ldq );
//	info = LAPACKE_dorgqr_work( LAPACK_COL_MAJOR, ml, n, n, Qii, ldq, tau, work, lwork );


// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
// // // // //   	       		  This block computes the Cholesky QR					     //
// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //

	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml, n, Aii, lda, Qii, ldq ); 

	cblas_dsyrk( CblasColMajor, CblasUpper, CblasTrans, n, ml, 1.0e+00, Qii, ldq, 0e+00, R, n );
	info = LAPACKE_dpotrf( LAPACK_COL_MAJOR, 'U', n, R, n ); 
	cblas_dtrsm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, ml, n, 1.0e+00, R, n, Qii, ldq );

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//				This block computes the LU-Factorization and changes sign in R & Q 		      //
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml, n, Qii, ldq, Aii, lda );

	for( i1 = 0; i1 < n ; i1++){
		if ( fabs(1 - Aii[ i1 + i1*lda ]) < fabs(1 + Aii[ i1 + i1*lda ]) ){

			for( j1 = i1; j1 < n ; j1++) R[ i1 + j1*n ] = - R[ i1 + j1*n ];
			for( j1 = 0; j1 < ml ; j1++) Qii[ j1 + i1*ldq ] = - Qii[ j1 + i1*ldq ];

		} else {

			for( j1 = 0; j1 < ml ; j1++) Aii[ j1 + i1*lda ] = - Aii[ j1 + i1*lda ];			

		}

		Aii[ i1 + i1*lda ] += 1.0e+00; 
		for( j1 = i1+1; j1 < ml; j1++ ) Aii[ j1 + i1*lda ] = Aii[ j1 + i1*lda ] / Aii[ i1 + i1*lda ];
		cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, ml-i1-1, n-i1-1, 1, (-1.0e+00), Aii + (i1+1) + i1*lda, lda, Aii + i1 + (i1+1)*lda, lda, (1.0e+00), Aii + (i1+1) + (i1+1)*lda, lda ); 

	}

// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
// // // // //   Trying to put back the good T now
// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
/*

	double *TTTki;
	int vb, k, j, itlo;

	k = i % mt;
	itlo = i % mt; // if( itlo == 0 ) itlo = mt;

	vb = mt - k; if ( vb > n ) vb = n;
	printf("vb = %d,\n", vb);

	TTTki = TTT + itlo + i*llldddttt;

	cblas_dtrsm( CblasColMajor, CblasRight, CblasLower, CblasTrans, CblasUnit, vb, vb, 1.0e+00, Asave, ml, TTTki, llldddttt );

	TTTki = TTTki + itlo + vb*llldddttt;
	j = i + vb;
	itlo = ((itlo + vb) % mt);

	if( j + mt >= i + n ) vb = n - ( j - i ); else vb = mt;

	while( vb != 0 ){

		cblas_dtrsm( CblasColMajor, CblasRight, CblasLower, CblasTrans, CblasUnit, vb, vb, 1.0e+00, Asave, ml, TTTki, llldddttt );

		itlo = ((itlo + vb) % mt);		
		TTTki = TTTki + itlo + vb*llldddttt;
		j += vb;
		if( j + mt >= i + n ) vb = n - ( j - i ); else vb = mt;
	}

*/

// // // // //
// // // // //     Put back the good R
// // // // //

	info = lila_dlarft_w03( m, n, i, mt, A, lda, TTT, llldddttt, tau);


	for( i1 = 0; i1 < ml ; i1++){
		for( j1 = 0; j1 < n ; j1++){
			if( i1 <= j1 ) Aii[ i1 + j1*lda ] = R[ i1 + j1*n ];
		}
	}
	
	free( R );
	free( tau );
//	free( Asave );
//	free( Qsave );

	return 0;

}
