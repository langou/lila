#include "lila.h"

int lila_dorghr( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork ){

	int info; 
	double *Aii, *Qii;
	int ml;
	double *Tki;
	int vb, k, j;
	int *S;
	int i1, j1;

	S = (int *) malloc(mt * sizeof(int));

	ml = m - i;
	k = i % mt;

	Aii = A + i*lda + i;
	Qii = Q + i*ldq + i;
	Tki = T + k + i*ldt;
	
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'L', ml-1, n, Qii+1, ldq, Aii+1, lda );

	vb = mt - k; if ( vb > n ) vb = n;

	lila_dorgh2( ml, vb, mt, Aii, lda, Tki, ldt, Qii, ldq, work, lwork, S );


	for( i1 = 0; i1 < n ; i1++){

		if ( S[ i1 ] == -1 ){

			for( j1 = i1; j1 < n ; j1++) Aii[ i1 + j1*lda ] = - Aii[ i1 + j1*lda ];
			for( j1 = 0; j1 < ml ; j1++) Qii[ j1 + i1*ldq ] = - Qii[ j1 + i1*ldq ];
		}

	}

	j = i + vb;

double *Qij, *Qjj, *Ajj, *Aji, *workjj;
	
	Tki = T + j*ldt;
	Qij = Q + j*lda + i;
	Ajj = Aii + vb*(1+lda);
	Aji = A + i*lda + j;

	Qii += vb*(1+ldq);

	Qij = Q + i + j*ldq;
	Qjj = Q + j + j*ldq;

	if( j+mt >= i+n ) vb = n-(j-i); else vb = mt;

	while( vb != 0 ){

		workjj = work + j*vb;

		info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', j, vb, Qij, ldq, work, j ); 

		info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', vb, vb, Qjj, ldq, workjj, vb ); 

		cblas_dtrsm( CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, j, vb, 1.0e+00, Aii, lda, work, j );

		cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, vb, vb, j, (-1.0e+00), Aji, lda, work, j, (1.0e+00), workjj, lda );


		lila_dorgh2( ml, vb, mt, Ajj, lda, Tki, ldt, Qii, ldq, work, lwork, S );


	for( i1 = 0; i1 < n ; i1++){

		if ( S[ i1 ] == -1 ){

			for( j1 = i1; j1 < n ; j1++) Aii[ i1 + j1*lda ] = - Aii[ i1 + j1*lda ];
			for( j1 = 0; j1 < ml ; j1++) Qii[ j1 + i1*ldq ] = - Qii[ j1 + i1*ldq ];
		}

	}


		j += vb;

		Tki += ( vb*ldt );
		Ajj += vb*(1+lda);
		Qii += vb*(1+ldq);

		if( j+mt >= i+n ) vb = n-(j-i); else vb = mt;
	}

	free(S);

	return 0;

}
