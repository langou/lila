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

	printf("\n entering now \n");

	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'L', ml-1, vb, Qii+1, ldq, Aii+1, lda ); 
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'U', vb, vb, Qii, ldq, Tki, ldt );
 
	j = i;
//	lila_dorgh2( m, vb, 0, mt, Aii, lda, Tki, ldt, NULL, -1, work, lwork, S );
	lila_dorgh2( m, n, j, mt, A, lda, T, ldt, Q, ldq, work, lwork, S );

	printf("\n");

	free( zork );
	free(TTT);
	free(Asave);

	return 0;

}
