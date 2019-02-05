#include "lila.h"

int lila_dorghr( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork, int *S ){

	int info; 
	double *Aii, *Qii;
	int ml;
	double *Tki;
	int vb, k, j;
	int i1, j1;
	double *TTT, *Asave;
	double *zork;

	int *Svb;
	zork = (double *) malloc(n * n * sizeof(double));
	TTT = (double *) malloc(n * n * sizeof(double));
	Asave = (double *) malloc(lda * n * sizeof(double));

	ml = m - i;
	k = i % mt;
	vb = mt - k; if ( vb > n ) vb = n;

	Aii = A + i*lda + i;
	Qii = Q + i*ldq + i;
	Tki = T + k + i*ldt;

//	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, n, Aii, lda, Asave + i + i*lda, lda ); 
//	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'L', ml-1, n, Qii+1, ldq, Asave + i+1 + i*lda, lda ); 
//	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'U', n, n, Qii, ldq, TTT, n );
//	lila_dorgh2( m, n, 0, mt, Asave + i + i*lda, lda, TTT, n, NULL, -1, work, lwork, S );
	printf("\n entering now \n");

	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'L', ml-1, vb, Qii+1, ldq, Aii+1, lda ); 
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'U', vb, vb, Qii, ldq, Tki, ldt );
 
	j = i;
//	lila_dorgh2( m, vb, 0, mt, Aii, lda, Tki, ldt, NULL, -1, work, lwork, S );
	lila_dorgh2( m, n, j, mt, A, lda, T, ldt, Q, ldq, work, lwork, S );

	printf("\n");






/*


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

		info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', j-i, vb, Qij, ldq, work, j-i ); 

		info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m-j, vb, Qjj, ldq, zork, m-j ); 

		// the -1 is a cheat
		cblas_dtrsm( CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, j-i, vb, -1.0e+00, Aii, lda, work, j-i ); 

		cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, m-j, vb, j-i, (-1.0e+00), Aji, lda, work, j-i, (+1.0e+00), zork, m-j ); 

		info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'L', m-j-1, vb, zork+1, m-j, Ajj+1, lda ); 

		info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'U', vb, vb, Qjj, ldq, T0j, ldt ); 

		lila_dorgh2( m-j, vb, 0, mt, Ajj, lda, T0j, ldt, NULL, -1, work, lwork, Svb );

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
*/
	free( zork );
	free(TTT);
	free(Asave);

	return 0;

}