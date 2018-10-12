#include "lila.h"

int lila_dge_qr_ormqrf_w03( int m, int n, int k, int i, int j, int mt, double *A, int lda, double *T, int ldt, double *work, int lwork ){

	double *Aii, *Tii, *Aij;
	int ml, vb;
	int ldwork;
	int iii, jjj; 

	int jalo;
	int not_done, jj;

	vb= mt - (i % mt );
	if ( vb > k ) vb = k;

	jalo = i;
	
	Aii = A + i + i*lda;
	Aij = A + i + j*lda;
	Tii = T + i + i*ldt;

	ldwork = mt;

	ml = m - i;

	not_done = 1;

	jj = 1;

	while( not_done == 1 ){

		printf("i=%d, mt=%d   (i %% mt)= %d, vb = mt - (i %% mt ) = %d, jalo = %d, j = %d\n",i,mt,(i % mt),vb,jalo,j);

		for( jjj = 0; jjj < n; jjj++ ){
			for( iii = 0; iii < vb; iii++ ){
				work[ iii + jjj * ldwork ] = Aij[ iii + jjj * lda  ];
			}
		}

 		cblas_dtrmm( CblasColMajor, CblasLeft, CblasLower, CblasTrans, CblasUnit, vb, n, (1.0e+00), Aii, lda, work, ldwork );
 		cblas_dgemm( CblasColMajor, CblasTrans, CblasNoTrans, vb, n, ml-vb, (1.0e+00), Aii+vb, lda, Aij+vb, lda, (1.0e+00), work, ldwork );
		cblas_dtrmm( CblasColMajor, CblasLeft, CblasUpper, CblasTrans, CblasNonUnit, vb, n, (1.0e+00), Tii, ldt, work, ldwork );
		cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, ml-vb, n, vb, (-1.0e+00), Aii+vb, lda, work, ldwork, (1.0e+00), Aij+vb, lda );
		cblas_dtrmm( CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, vb, n, (+1.0e+00), Aii, lda, work, ldwork );

		for( jjj = 0; jjj < n; jjj++ ){
			for( iii = 0; iii < vb; iii++ ){
				Aij[ iii + jjj * lda  ] -= work[ iii + jjj * ldwork ];
			}
		}

		if ( jj + vb - 1 == k ) {

			not_done = 0; 

		} else if ( jj + vb - 1 > k ) {

			printf("defensive programming, we should never have been there, abort\n"); return 0;

		} else {

			ml = ml-vb;

			jalo = jalo+vb;

			jj = jj + vb;
 
			if ( ( jj + mt - 1 ) <= k ) vb = mt; else vb = k-jj+1;

			Aij = A + jalo + jalo*lda;
			Aij = A + jalo + j*lda;
			Tii = T + jalo + jalo*ldt;

		}

	}

	return 0;

}
