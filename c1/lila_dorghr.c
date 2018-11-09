#include "lila.h"

int lila_dorghr( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork ){

	int info; 
	double *Aii, *Qii;
	int ml;
	double *Tki;
	int vb, k, j;
	int *S;
	int i1, j1;
	double *TTT, *Asave;
	double *zork;

	int *Svb;
	zork = (double *) malloc(n * n * sizeof(double));
	S = (int *) malloc(n * sizeof(int));
	TTT = (double *) malloc(n * n * sizeof(double));
	Asave = (double *) malloc(lda * n * sizeof(double));


	printf("\n entering now \n");

	ml = m - i;
	k = i % mt;
	vb = mt - k; if ( vb > n ) vb = n;

	Aii = A + i*lda + i;
	Qii = Q + i*ldq + i;
	Tki = T + k + i*ldt;

	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, n, Aii, lda, Asave + i + i*lda, lda ); 

	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'L', ml-1, n, Qii+1, ldq, Asave + i+1 + i*lda, lda ); 
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'U', n, n, Qii, ldq, TTT, n ); // Copy Q into T for dorgh2
	lila_dorgh2( ml, n, mt, Asave + i + i*lda, lda, TTT, n, NULL, -1, work, lwork, S );

//	for( j1 = i; j1 < m; j1++){
//		for( i1 = i; i1 < i+vb; i1++){
//			if(i1<j1) printf(" %+9.4e ", Asave[ j1 + i1*lda ]); else printf("  0.00000    ");
//		}
//		printf("\n");
//	}
//	printf("\n");


	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'L', ml-1, vb, Qii+1, ldq, Aii+1, lda ); 
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'U', vb, vb, Qii, ldq, Tki, ldt ); // Copy Q into T for dorgh2

	lila_dorgh2( ml, vb, mt, Aii, lda, Tki, ldt, NULL, -1, work, lwork, S );

	printf("\n");

	for( i1 = 0; i1 < vb ; i1++){

		if ( S[ i1 ] == -1 ){

			//for( j1 = i1; j1 < n ; j1++) Aii[ i1 + j1*lda ] = - Aii[ i1 + j1*lda ];
			//for( j1 = 0; j1 < ml ; j1++) Qii[ j1 + i1*ldq ] = - Qii[ j1 + i1*ldq ];
		}

	}

	Svb = S + vb;

	j = i + vb;
	ml -= vb;

//////////////////////////////////////////////////

	double *Qij, *Qjj, *Ajj, *Aji, *T0j;

	T0j = T + j*ldt;
	Qij = Qii + vb*ldq;
	Qjj = Qii + vb*(1+ldq);
	Ajj = Aii + vb*(1+lda);
	Aji = Aii + vb;

	if( j+mt >= i+n ) vb = n-(j-i); else vb = mt;

	while( vb != 0 ){

//printf("\n     We're  in  the  while %d, j=%d, %d, m=%d, vb=%d  \n\n", i, j, k, m, vb);


		info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', j-i, vb, Qij, ldq, work, j-i ); // Copy top part of Qi into work

		info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m-j, vb, Qjj, ldq, zork, m-j ); // Copy lower part of Qj into zork

		// the -1 is a cheat
		cblas_dtrsm( CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, j-i, vb, -1.0e+00, Aii, lda, work, j-i ); // Update work with L \ Qi

		cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, m-j, vb, j-i, (-1.0e+00), Aji, lda, work, j-i, (+1.0e+00), zork, m-j ); // Update zork with Qj - L*Qi

		info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'L', m-j-1, vb, zork+1, m-j, Ajj+1, lda );  // Copy the updated part of zork into the lower part of A

		info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'U', vb, vb, Qjj, ldq, T0j, ldt ); // Copy Q into T for dorgh2

		lila_dorgh2( m-j, vb, mt, Ajj, lda, T0j, ldt, NULL, -1, work, lwork, Svb );


	//for( i1 = 0; i1 < vb ; i1++){

		//if ( Svb[ i1 ] == -1 ){

			//for( j1 = i1; j1 < n-j ; j1++) Ajj[ i1 + j1*lda ] = - Ajj[ i1 + j1*lda ];
			//for( j1 = 0; j1 < m-i ; j1++) Qij[ j1 + i1*ldq ] = - Qij[ j1 + i1*ldq ];

		//}

	//}

		j += vb;

		T0j += ( vb*ldt );
		Aji += vb;
		Ajj += vb*(1+lda);
		Qij += vb*ldq;
		Qjj += vb*(1+ldq);

		Svb += vb;

		ml -= vb;

		if( j+mt >= i+n ) vb = n-(j-i); else vb = mt;
	
	}

	free( zork );
	free(S);
	free(TTT);
	free(Asave);

	return 0;

}
