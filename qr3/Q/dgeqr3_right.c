#include "qr3.h"

int dgeqr3_right( int m, int n, double *A, int lda, double *T, int ldt ){

	/* Only the upper part of T is referenced */

	int n1, n2, i,j;
	double *W;
	int ldw, info;

	if ( n == 1){
		info = LAPACKE_dlarfg_work( m, (&A[0]), &(A[1]), 1, &(T[0]));
	} else {

	n1 = n/2;
	n2 = n-n1;

	dgeqr3_left( m, n1, &(A[0]), lda, T, ldt );

		/* ------------------------------------------------------------------------------------- */
		/*     routine  | computation                     | number of floating-point operations  */
		/* ------------------------------------------------------------------------------------- */
		/* (1) DTRMM    | W = Y1^T * C1                   | k * ( k - 1 ) * n                    */
		/* (2) DGEMM    | W = W + Y2^T * C2               | 2 * k * ( m - k ) * n                */
		/* (3) DTRMM    | W = T^T * W                     | k^2 * n                              */
		/* (4) DGEMM    | W = C2 - Y2 * W                 | 2 * k * ( m - k ) * n                */
		/* (5) DTRMM    | W = Y1 * W                      | k * ( k - 1 ) * n                    */
		/* (6) matrix - | C1 = C1 - W                     | k * n                                */
		/* ------------------------------------------------------------------------------------- */
		/*     total    | W = - T1 * ( Y1^T * Y2 ) * T2   | k * n * ( 4 * m - k - 1 )            */
		/* ------------------------------------------------------------------------------------- */

		ldw=ldt;
		W = &(T[n1*ldt]);
		//W = (double *)malloc(ldw*n*sizeof(double)) ; 
		for (i=0;i<n1;i++) for (j=0;j<n2;j++) W[i+j*ldw] = A[(j+n1)*lda+i];
		cblas_dtrmm ( CblasColMajor, CblasLeft, CblasLower, CblasTrans, CblasUnit, n1, n2, 1.0e+00, &(A[0]), lda, W, ldw); 
		/*      (2) DGEMM    | W = W + Y2^T * C2               | 2 * k * ( m - k ) * n           */
		cblas_dgemm( CblasColMajor, CblasTrans, CblasNoTrans, n1, n2, m-n1, +1.0e+00, &(A[n1]), lda, &(A[n1+n1*lda]), lda, +1.0e+00, W, ldw);
		/*      (3) DTRMM    | W = T^T * W                     | k^2 * n                         */
		cblas_dtrmm ( CblasColMajor, CblasLeft, CblasUpper, CblasTrans, CblasNonUnit, n1, n2, 1.0e+00, &(T[0]), ldt, W, ldw); 
		/*      (4) DGEMM    | C2 = C2 - Y2 * W                | 2 * k * ( m - k ) * n           */
		cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, m-n1, n2, n1, -1.0e+00, &(A[n1]), lda, W, ldw, +1.0e+00, &(A[n1+n1*lda]), lda);
		/*      (5) DTRMM    | W = Y1 * W                      | k * ( k - 1 ) * n               */
		cblas_dtrmm ( CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, n1, n2, 1.0e+00, &(A[0]), lda, W, ldw); 
		/*      (6) matrix - | C1 = C1 - W                     | k * n                           */
		for (j=0;j<n2;j++) for (i=0;i<n1;i++) A[(j+n1)*lda+i] -= W[j*ldw+i];
		//free(W);

	dgeqr3_right( m-n1, n2, &(A[n1*lda+n1]), lda, &(T[n1+n1*ldt]), ldt );


	}

	return 0;

}