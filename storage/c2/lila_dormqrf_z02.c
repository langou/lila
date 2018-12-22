#include "lila.h"

int lila_dormqrf_z02( int m, int n, int k, int i, int j, int mt, double *A, int lda, double *B, int ldb, double *T, int ldt, double *work, int lwork ){

	double *Aii, *Tii, *Bij;
	int ml;
	int ldwork;
	int ii, jj; 

	Aii = A + i + i*lda;
	Bij = B + i + j*ldb;
	Tii = T + i + i*ldt;
	ml = m - i;
	ldwork = k;

	for( jj = 0; jj < n; jj++ ){
	for( ii = 0; ii < k; ii++ ){
	work[ ii + jj * ldwork ] = Bij[  ii + jj * ldb  ];
	}}

//	printf("ORMQRF  we need work to be of size %2dx%2d\n",k,n);

 	cblas_dtrmm( CblasColMajor, CblasLeft, CblasLower, CblasTrans, CblasUnit, k, n, (1.0e+00), Aii, lda, work, ldwork );
 	cblas_dgemm( CblasColMajor, CblasTrans, CblasNoTrans, k, n, m-k-i, (1.0e+00), Aii+k, lda, Bij+k, ldb, (1.0e+00), work, ldwork );
	cblas_dtrmm( CblasColMajor, CblasLeft, CblasUpper, CblasTrans, CblasNonUnit, k, n, (1.0e+00), Tii, ldt, work, ldwork );
	cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, m-k-i, n, k, (-1.0e+00), Aii+k, lda, work, ldwork, (1.0e+00), Bij+k, ldb );
	cblas_dtrmm( CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, k, n, (1.0e+00), Aii, lda, work, ldwork );

	for( jj = 0; jj < n; jj++ ){
	for( ii = 0; ii < k; ii++ ){
	Bij[  ii + jj * ldb  ] -= work[ ii + jj * ldwork ];
	}}

	return 0;

}
