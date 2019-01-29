#include "lila.h"

int lila_dgeqrf_w03_mt_v02( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork ){

	double *Tki, *Aii, *Qii;
	int ml, vb, j, k, info;
	int *S;
	
	S  = (int *) malloc(n * sizeof(int));
	
	j  = i;
	k  = i % mt;
	ml = m - i;
	vb = mt - ( i % mt); if ( vb > n ) vb = n;

	Aii = A + i        + i*lda;
	Qii = Q + i        + i*ldq;
	Tki = T + (i % mt) + i*ldt;

	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml, n, Aii, lda, Qii, ldq ); 
	cblas_dsyrk( CblasColMajor, CblasUpper, CblasTrans, n, ml, 1.0e+00, Qii, ldq, 0e+00, Aii, lda );
	info = LAPACKE_dpotrf( LAPACK_COL_MAJOR, 'U', n, Aii, lda ); 
	cblas_dtrsm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, ml, n, 1.0e+00, Aii, lda, Qii, ldq );

	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml-vb, vb, Qii+vb, ldq, Aii+vb, lda ); 
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', vb, vb, Qii, ldq, Tki, ldt ); 
	lila_dorgh2_3( ml, vb, Aii, lda, Tki, ldt, Qii, ldq, work, lwork, S );

	j  += vb;
	ml -= vb;

	Aii += vb*(1+lda);
	Qii += vb*(1+ldq);
	Tki += vb;

	if( j + mt >= i + n ) vb = n - ( j - i ); else vb = mt;

	while( vb != 0 ){

		info = lila_ormhr2_w03_hr( m, vb, i, j, 0, mt, A, lda, T, ldt, Q, ldq, work, lwork, S );
		info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml-vb, vb, Qii+vb, ldq, Aii+vb, lda ); 
		info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', vb, vb, Qii, ldq, Tki, ldt ); 
		lila_dorgh2_3( ml, vb, Aii, lda, Tki, ldt, Qii, ldq, work, lwork, S );

		Aii += vb*(1+lda);
		Qii += vb*(1+ldq);
		Tki += vb;
		j   += vb;
		ml  -= vb;

		if( j + mt >= i + n ) vb = n - ( j - i ); else vb = mt;

	}

	free(S);

	return 0;

}
