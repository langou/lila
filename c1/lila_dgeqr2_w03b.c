#include "lila.h"

int lila_dgeqr2_w03b( int m, int n, int i, int mt, double *A, int lda, double *TTT, int llldddttt, double *Q, int ldq, double *work, int lwork ){

	int info ; 
	double *tau=NULL;
	double *Aii, *Qii, *TTTii;
	int ml;

	tau = (double *) malloc( n * sizeof(double));

	Aii = A + i*lda + i;
	Qii = Q + i*ldq + i;
	TTTii = TTT + i*llldddttt + i;
	
	printf("\n");
	printf("m = %d, i = %d, n = %d,",m,i,n);
	printf("\n");
	printf("\n");
	ml = m - i;

	double *Asave, *Qsave, *R, *Asaveii;
	int ldasave;
	int i1, j1;

	Asave = (double *) malloc( ml * n * sizeof(double));
	Qsave = (double *) malloc( ml * n * sizeof(double));
	R = (double *) malloc( n * n * sizeof(double));



// // // // //
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml, n, Aii, lda, Qsave, ml ); 
// // // // //


// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
// // // // //    This block is to make sure we still have a working code. I copy everything else to work with and will place pieces back to see if I can get it working.
// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //


	info = LAPACKE_dgeqrf_work( LAPACK_COL_MAJOR, ml, n, Aii, lda, tau, work, lwork );
	info = lila_dlarft_w03( m, n, i, mt, A, lda, TTT, llldddttt, tau);
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml, n, Aii, lda, Qii, ldq );
	info = LAPACKE_dorgqr_work( LAPACK_COL_MAJOR, ml, n, n, Qii, ldq, tau, work, lwork );

	printf("\n");   


// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
// // // // //   	       		  This block computes the Cholesky QR					     //
// // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //


	cblas_dsyrk( CblasColMajor, CblasUpper, CblasTrans, n, ml, 1.0e+00, Qsave, ml, 0e+00, R, n );
	info = LAPACKE_dpotrf( LAPACK_COL_MAJOR, 'U', n, R, n ); 
	cblas_dtrsm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, ml, n, 1.0e+00, R, n, Qsave, ml );


// // // // //
// // // // //
// // // // //

	printf("R First = \n");
	for( i1 = 0; i1 < n ; i1++){
	for( j1 = 0; j1 < n ; j1++){
		printf("% 5.3e, ",R[ i1 + j1*n]);
	}
	printf("\n");
	}
	printf("\n");

// 			This one was the cheat Julien started this script using
//
//	for( i1 = 0; i1 < n ; i1++){
//		if ( R[ i1 + i1*n ] * Aii[ i1 + i1*lda ] < 0 ){
//			for( j1 = 0; j1 < n ; j1++) R[ i1 + j1*n ] = - R[ i1 + j1*n ];
//			for( j1 = 0; j1 < ml ; j1++) Qii[ j1 + i1*ldq ] = - Qii[ j1 + i1*ldq ];
//		}
//	}


	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml, n, Qsave, ml, Asave, ml );

	for( i1 = 0; i1 < n ; i1++){
		if ( fabs(1 - Asave[ i1 + i1*ml ]) < fabs(1 + Asave[ i1 + i1*ml ]) ){
			for( j1 = i1; j1 < n ; j1++) R[ i1 + j1*n ] = - R[ i1 + j1*n ];
			for( j1 = 0; j1 < ml ; j1++) Qsave[ j1 + i1*ml ] = - Qsave[ j1 + i1*ml ];
		} else {
			for( j1 = 0; j1 < ml ; j1++) Asave[ j1 + i1*ml ] = - Asave[ j1 + i1*ml ];			
		}
		// This is the LU putting it right into Asave, L is lower unit and U is triu(Asave)
		Asave[ i1 + i1*ml ] += 1; 
		for( j1 = i1; j1 < ml; j1++ ) Asave[ j1 + i1*ml ] = Asave[ j1 + i1*ml ] / Asave[ i1 + i1*ml ]; 
		cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, ml-i1, n-i1-1, 1, (-1.0e+00), Asave + i1 + (i1-1)*ml, ml, Asave + (i1-1) + i1*ml, ml, (1.0e+00), Asave + i1 + i1*ml, ml ); 
	}

	printf("R Second (for sign changes) = \n");
	for( i1 = 0; i1 < n ; i1++){
	for( j1 = 0; j1 < n ; j1++){
		printf("% 5.3e, ", R[ i1 + j1*n ]);
	}
	printf("\n");
	}

	printf("\nSee if A is the same as R \n");
	for( i1 = 0; i1 < n; i1++){
	for( j1 = 0; j1 < n; j1++){
		if( i1 <= j1 ) printf("% 5.3e, ", Aii[ i1 + j1*lda ]); else printf("% 5.3e, ", 0.00e00);
	}
	printf("\n");   
	}
	// After comparing R and triu(A) we can see that they are the same but aren't capturing the good sign changes.


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



/////   The below commented section puts zeros in TTT, copies R in TTT, and compares with the original T (first print). I realized it didn't follow the mt 
////     structure bby comparing the original T with TTT so I started working on the above block to get the good T.
/*
// set T to zeros, copy R in T, trsm with Lower, Transpose, Unit, Right

	printf("\nThe good TTT (before putting zeros) = \n");   
	for( i1 = 0; i1 < vb ; i1++){
	for( j1 = 0; j1 < n ; j1++){
		printf("% 5.3e, ", TTTki[ i1 + j1*llldddttt ] );
		}
	printf("\n");   
	}
	printf("\n");   

	for( i1 = 0; i1 < n ; i1++){
	for( j1 = 0; j1 < n ; j1++){
		TTTii[ i1 + j1*llldddttt ]= 0.0e+00;
		}
	}

	printf("\nTTT (with new R input within) = \n");   
	for( i1 = 0; i1 < n ; i1++){
	for( j1 = 0; j1 < n ; j1++){
		TTTii[ i1 + j1*llldddttt ]= R[ i1 + j1*n];
		printf("% 5.3e, ", TTTii[ i1 + j1*llldddttt ] );
		}
	printf("\n");   
	}

//	cblas_dtrsm( CblasColMajor, CblasRight, CblasLower, CblasTrans, CblasUnit, ml, n, 1.0e+00, Asave, ml, TTTii, llldddttt );

*/
// // // // //
// // // // //     Put back the good R
// // // // //

	for( i1 = 0; i1 < n ; i1++){
	for( j1 = i1; j1 < n ; j1++){
		if( i1 <= j1 ) Asave[ i1 + j1*ml ] = R[ i1 + j1*n ];  
	}
	}

	free( R );
	free( tau );
	free( Asave );
	free( Qsave );

	return 0;

}
