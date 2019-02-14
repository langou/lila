#include "lila.h"

int lila_dgeqrf_w03_mt_hr( int panel, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork ){

	double *Tki, *Aii, *Qii, *S;
	int ml, vb, j, info;
	
//	S  = (int *) malloc((n+i) * sizeof(int));
	S = work;	

	j  = i;
	ml = m - i;
	vb = mt - ( i%mt ); if ( vb > n ) vb = n;

	Aii = A + i        + i*lda;
	Qii = Q + i        + i*ldq;
	Tki = T + ( i%mt ) + i*ldt;

	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml, n, Aii, lda, Qii, ldq ); 
	cblas_dsyrk( CblasColMajor, CblasUpper, CblasTrans, n, ml, (+1.0e+00), Qii, ldq, (+0e+00), Aii, lda );
	info = LAPACKE_dpotrf( LAPACK_COL_MAJOR, 'U', n, Aii, lda ); 
	cblas_dtrsm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, ml, n, (+1.0e+00), Aii, lda, Qii, ldq );

	while( vb != 0 ){

		info = lila_ormhr2_w03_hr( m, vb, i, j, -1, mt, A, lda, T, ldt, Q, ldq, work, lwork, S );
		j += vb;
		if( j + mt >= i + n ) vb = n - ( j - i ); else vb = mt;

	}

//	free(S);

	return 0;

}
