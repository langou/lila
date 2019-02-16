#include "lila.h"

int lila_wsq_dormqrbz_w03( int m, int n, int k, int i, int j, int mt, double *A, int lda, double *Q, int ldq, double *T, int ldt, double *work, int lwork ){

	int vb, jj, lwork1, not_done, ldwork, ml;

	vb = ((i+k) % mt); if (vb == 0) vb = mt; if ( vb > k ) vb = k;
	lwork1 = 0;
//	ldwork = k;
	
	if( lwork1 < n*ldwork ) lwork1 = n*ldwork;
//	if( lwork1 < k ) lwork1 = k;

	ml = m - (i+k-vb);
	not_done = 1;
	jj = 1;

//	cblas_dgemm( CblasColMajor, CblasTrans, CblasNoTrans, vb, n, ml-vb, (1.0e+00), Aii+vb, lda, Qij+vb, ldq, (0.0e+00), work, ldwork );
//	cblas_dtrmm( CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, vb, n, (1.0e+00), Tii, ldt, work, ldwork );
//	LAPACKE_dlacpy( LAPACK_COL_MAJOR, 'A', vb, n, work, ldwork, Qij, ldq );
//	cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, ml-vb, n, vb, (-1.0e+00), Aii+vb, lda, work, ldwork, (1.0e+00), Qij+vb, ldq );

	while ( not_done == 1 ){
		if ( jj + vb - 1 == k ){ 
			not_done = 0; 
		} else if ( jj + vb - 1 > k ){
			printf("defensive programming, we should never have been there, abort\n");
		} else{

			if( lwork1 < n*vb ) lwork1 = n*vb;
			if( lwork1 < ml-vb ) lwork1 = ml-vb;

			jj += vb;
			if ( ( jj + mt - 1 ) <= k ) vb = mt; else { vb = k - jj + 1;  }
			ml += vb;
// 			cblas_dgemm( CblasColMajor, CblasTrans, CblasNoTrans, vb, n, ml-vb, (1.0e+00), Aii+vb, lda, Qij+vb, ldq, (0.0e+00), work, ldwork );
			// The inner dimension for Q*work is ml-vb. Does that mean this neead to be the size?

//			cblas_dtrmm( CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, vb, n, (1.0e+00), Tii, ldt, work, ldwork );
//			LAPACKE_dlacpy( LAPACK_COL_MAJOR, 'A', vb, n, work, ldwork, Qij, ldq );
//			cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, ml-vb, n, vb, (-1.0e+00), Aii+vb, lda, work, ldwork, (1.0e+00), Qij+vb, ldq );
		}
	}

	return lwork1;

}
