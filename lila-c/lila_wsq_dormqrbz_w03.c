#include "lila.h"

int lila_wsq_dormqrbz_w03( int m, int n, int k, int i, int j, int mt, double *A, int lda, double *Q, int ldq, double *T, int ldt, double *work, int lwork ){

	double *Aii, *Qij, *Tii;
	int vb, jj, not_done, ldwork, ml;

	int itlo;
	int ithi;
	int jtlo;

	vb = ((i+k) % mt); if (vb == 0) vb = mt; ithi = vb; if ( vb > k ) vb = k;
	itlo = ithi - vb;
	jtlo = (i+k-vb);

	Aii = A + jtlo + jtlo*lda; 
	Qij = Q + jtlo + j*ldq;
	Tii = T + itlo + jtlo*ldt;

	ldwork = k;
	ml = m - (i+k-vb);
	not_done = 1;
	jj = 1;

	cblas_dgemm( CblasColMajor, CblasTrans, CblasNoTrans, vb, n, ml-vb, (1.0e+00), Aii+vb, lda, Qij+vb, ldq, (0.0e+00), work, ldwork );
	cblas_dtrmm( CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, vb, n, (1.0e+00), Tii, ldt, work, ldwork );
	LAPACKE_dlacpy( LAPACK_COL_MAJOR, 'A', vb, n, work, ldwork, Qij, ldq );
	cblas_dtrmm( CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, vb, n, (-1.0e+00), Aii, lda, Qij, ldq );
	cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, ml-vb, n, vb, (-1.0e+00), Aii+vb, lda, work, ldwork, (1.0e+00), Qij+vb, ldq );

	while ( not_done == 1 ){

		if ( jj + vb - 1 == k ){ 

			not_done = 0; 

		} else if ( jj + vb - 1 > k ){

			printf("defensive programming, we should never have been there, abort\n");

		} else{

			jj += vb;

			if ( ( jj + mt - 1 ) <= k ) vb = mt; else { vb = k - jj + 1; Tii += ( i % mt ) ; }

			ml += vb;

			Aii -= vb + vb*lda;
			Qij -= vb;
			Tii -= vb*ldt;

 			cblas_dgemm( CblasColMajor, CblasTrans, CblasNoTrans, vb, n, ml-vb, (1.0e+00), Aii+vb, lda, Qij+vb, ldq, (0.0e+00), work, ldwork );
			cblas_dtrmm( CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, vb, n, (1.0e+00), Tii, ldt, work, ldwork );
			LAPACKE_dlacpy( LAPACK_COL_MAJOR, 'A', vb, n, work, ldwork, Qij, ldq );
			cblas_dtrmm( CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, vb, n, (-1.0e+00), Aii, lda, Qij, ldq );
			cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, ml-vb, n, vb, (-1.0e+00), Aii+vb, lda, work, ldwork, (1.0e+00), Qij+vb, ldq );


		}

	}

	return 0;

}
