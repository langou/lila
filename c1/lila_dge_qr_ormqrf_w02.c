#include "lila.h"

int lila_dge_qr_ormqrf_w02( int m, int n, int k, int i, int j, int mt, double *A, int lda, double *T, int ldt, double *work, int lwork ){

	double *Aii, *Tii, *Aij;
	int ml;
	int ldwork;
	int ii, jj; 

	Aii = A + i + i*lda;
	Aij = A + i + j*lda;
	Tii = T + i + i*ldt;
	ml = m - i;
	ldwork = k;

	for( jj = 0; jj < n; jj++ ){
	for( ii = 0; ii < k; ii++ ){
	work[ ii + jj * ldwork ] = Aij[  ii + jj * lda  ];
	}}

 	cblas_dtrmm( CblasColMajor, CblasLeft, CblasLower, CblasTrans, CblasUnit, k, n, (1.0e+00), Aii, lda, work, ldwork );
 	cblas_dgemm( CblasColMajor, CblasTrans, CblasNoTrans, k, n, m-k-i, (1.0e+00), Aii+k, lda, Aij+k, lda, (1.0e+00), work, ldwork );
	cblas_dtrmm( CblasColMajor, CblasLeft, CblasUpper, CblasTrans, CblasNonUnit, k, n, (1.0e+00), Tii, ldt, work, ldwork );
	cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, m-k-i, n, k, (-1.0e+00), Aii+k, lda, work, ldwork, (1.0e+00), Aij+k, lda );
	cblas_dtrmm( CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, k, n, (1.0e+00), Aii, lda, work, ldwork );

	for( jj = 0; jj < n; jj++ ){
	for( ii = 0; ii < k; ii++ ){
	Aij[  ii + jj * lda  ] -= work[ ii + jj * ldwork ];
	}}

	return 0;

}

/*
	int info;
	double *V;
	int vl;
	double *tau;
	double normV_square;

	tau = (double *) malloc( k * sizeof(double));

	for(jj = 0, vl=ml-1, V = Aii+1; jj < k; jj++,vl--,V+=(lda+1) ){
		normV_square = ( 1.0e+00 ) + cblas_ddot( vl, V, 1, V, 1 );
		tau[jj] = ( 2.0e+00 ) / normV_square ;
	}

	info = LAPACKE_dormqr_work( LAPACK_COL_MAJOR, 'L', 'T', ml, n, k, Aii, lda, tau, Aij, lda, work, lwork );

	free( tau );
	return 0;
*/

/*
	int info;
	double *V;
	int vl;
	double *tau;
	double normV_square;

	ldwork = n;

	tau = (double *) malloc( k * sizeof(double));
	for(jj = 0, vl=ml-1, V = Aii+1; jj < k; jj++,vl--,V+=(lda+1) ){
		normV_square = ( 1.0e+00 ) + cblas_ddot( vl, V, 1, V, 1 );
		tau[jj] = ( 2.0e+00 ) / normV_square ;
	}

	info = LAPACKE_dlarft_work ( LAPACK_COL_MAJOR, 'F', 'C', ml, k, Aii, lda, tau, Tii, ldt);

	info = LAPACKE_dlarfb_work ( LAPACK_COL_MAJOR, 'L', 'T', 'F', 'C', ml, n, k, Aii, lda, Tii, ldt, Aij, lda, work, ldwork );

	free( tau );
	return 0;
*/
