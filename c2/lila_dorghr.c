#include "lila.h"

int lila_dorghr( int m, int n, int i, int j, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork, int *S ){

	double *Aj0, *Ajj, *Tij, *Ti0;
	int k, i1;


	printf("\n entering now \n");
	printf("i = %d, j = %d, n = %d, m = %d\n",i,j,n,m);

	i = 0;
	Aj0 = A + j;
	Ajj = A + j + j*lda;
	Tij = T + i + j*ldt;
	Ti0 = T + i;

	lila_ormhr_w0b( m, n, i, j, A, lda, T, ldt, Q, ldq, S );

	lila_dorgh2( m, n, i, j, mt, A, lda, T, ldt, Q, ldq, work, lwork, S );

//	Make-shift connect
	for( k = i; k < j-1; k++){
		for( i1 = j; i1 < j+n-1; i1++){
		//cblas_dgemm( CblasColMajor, CblasNoTrans, CblasTrans, j-0, n, j, (-1.0e+00), Ti0, ldt, Aj0, lda, (+1.0e+00), Tij, ldt ); 
		//cblas_dtrsm( CblasColMajor, CblasRight, CblasLower, CblasTrans, CblasUnit, j-0, n, (+1.0e+00), Ajj, lda, Tij, ldt ); 
		}
	}

//	printf("\n");
//    for( k = 0; k < j+n; k++){
//	for( i1 = 0; i1 < j+n; i1++){
//	   printf(" %+5.2f ", T[ k + i1*ldt ] );
//	} 
//	printf("\n");
//   } 

	return 0;

}
