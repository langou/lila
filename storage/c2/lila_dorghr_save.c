#include "lila.h"

int lila_dorghr( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork ){

	int info ; 
	double *Aii, *Ajj, *Qii, *R;
	int ml, i1, j1;

	R = (double *) malloc(n * n * sizeof(double));

	Aii = A + i*lda + i;
	Qii = Q + i*ldq + i;
	
	ml = m - i;

	for( i1 = 0; i1 < ml ; i1++){
		for( j1 = 0; j1 < n ; j1++){
			if( i1 <= j1 ) R[ i1 + j1*n ] = Aii[ i1 + j1*lda ];
		}
	}

	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml, n, Qii, ldq, Aii, lda );

//	lila_dorgh2( m, n, i, mt, A, lda, Q, ldq, work, lwork );

	double Aii_val;

	for( i1 = 0; i1 < n ; i1++){
		if ( fabs(1 - Aii[ i1 + i1*lda ]) < fabs(1 + Aii[ i1 + i1*lda ]) ){

			for( j1 = i1; j1 < n ; j1++) R[ i1 + j1*n ] = - R[ i1 + j1*n ];
			for( j1 = 0; j1 < ml ; j1++) Qii[ j1 + i1*ldq ] = - Qii[ j1 + i1*ldq ];

		} else {

			for( j1 = 0; j1 < ml ; j1++) Aii[ j1 + i1*lda ] = - Aii[ j1 + i1*lda ];			

		}

		Aii_val = Aii[ i1 + i1*lda ] + 1.0e+00; 
		Aii[ i1 + i1*lda ] = Aii_val; 

		for( j1 = i1+1; j1 < ml; j1++ ) Aii[ j1 + i1*lda ] = Aii[ j1 + i1*lda ] / Aii[ i1 + i1*lda ];
		cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, ml-i1-1, n-i1-1, 1, (-1.0e+00), Aii + (i1+1) + i1*lda, lda, Aii + i1 + (i1+1)*lda, lda, (1.0e+00), Aii + (i1+1) + (i1+1)*lda, lda ); 

	}

//	info = lila_dlarft_w03_b( m, n, i, mt, A, lda, T, ldt );

	double *Tki;
	int vb, k, j;

	k = i % mt;

	Aii = A + i + i*lda;
	Tki = T + k + i*ldt;

	vb = mt - k; if ( vb > n ) vb = n;

	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'U', vb, vb, Aii, lda, Tki, ldt ); 
	cblas_dtrsm( CblasColMajor, CblasRight, CblasLower, CblasTrans, CblasUnit, vb, vb, 1.0e+00, Aii, lda, Tki, ldt );

	j = i + vb;
	
	Aii += vb*(1+lda);
	Tki = T + j*ldt;

	if( j+mt >= i+n ) vb = n-(j-i); else vb = mt;

	while( vb != 0 ){

		info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'U', vb, vb, Aii, lda, Tki, ldt ); 
		cblas_dtrsm( CblasColMajor, CblasRight, CblasLower, CblasTrans, CblasUnit, vb, vb, 1.0e+00, Aii, lda, Tki, ldt );

		j += vb;

		Tki += ( vb*ldt );
		Aii += vb*(1+lda);

		if( j+mt >= i+n ) vb = n-(j-i); else vb = mt;
	}

	Aii = A + i + i*lda;

	for( i1 = 0; i1 < ml ; i1++){
		for( j1 = 0; j1 < n ; j1++){
			if( i1 <= j1 ) Aii[ i1 + j1*lda ] = R[ i1 + j1*n ];
		}
	}

	free( R );

	return 0;

}
