#include "lila.h"

int lila_dgeqrf_w03_mt_v02( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork ){

	double *Tki, *Aii, *Qii, *AAs, *QQs;
	int ml, vb, j, info;
	int *S;
	
	S  = (int *) malloc((n+i) * sizeof(int));
	AAs = (double *) malloc(lda * (n+i) * sizeof(double));
	QQs = (double *) malloc(ldq * (n+i) * sizeof(double));
	
	j  = i;
	ml = m - i;
	vb = mt - ( i % mt); if ( vb > n ) vb = n;

	Aii = A + i        + i*lda;
	Qii = Q + i        + i*ldq;
	Tki = T + (i % mt) + i*ldt;

	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml, n, Aii, lda, Qii, ldq ); 
	cblas_dsyrk( CblasColMajor, CblasUpper, CblasTrans, n, ml, (+1.0e+00), Qii, ldq, (+0e+00), Aii, lda );
	info = LAPACKE_dpotrf( LAPACK_COL_MAJOR, 'U', n, Aii, lda ); 
	cblas_dtrsm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, ml, n, (+1.0e+00), Aii, lda, Qii, ldq );

// cheating /////////
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'U', n, n, Aii, lda, AAs, lda ); 
	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, n+i, Q, ldq, QQs, ldq ); 
/////////////////////

//	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', ml-vb, vb, Qii+vb, ldq, Aii+vb, lda ); 
//	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', vb, vb, Qii, ldq, Tki, ldt ); 
//	lila_dorgh2_3( ml, vb, Aii, lda, Tki, ldt, Qii, ldq, work, lwork, S );
	//j  += vb;
	//ml -= vb;
	//Aii = A + j        + j*lda;
	//Qii = Q + j        + j*ldq;
	//Tki = T + (j % mt) + j*ldt;
	//if( j + mt >= i + n ) vb = n - ( j - i ); else vb = mt;

	while( vb != 0 ){

		info = lila_ormhr2_w03_hr( m, vb, i, j, -1, mt, A, lda, T, ldt, Q, ldq, work, lwork, S );

		j += vb;
		if( j + mt >= i + n ) vb = n - ( j - i ); else vb = mt;

	}

//	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'U', n, n, AAs, lda, Aii, lda ); 
//	info = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, n+i, QQs, ldq, Q, ldq ); 

	free(S);
	free(AAs);
	free(QQs);

	return 0;

}
