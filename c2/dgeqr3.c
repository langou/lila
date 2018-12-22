#include "lila.h"

int dgeqr3( int m, int n, double *A, int lda, double *T, int ldt ){

	/* Only the upper part of T is referenced */

	int n1, n2, i,j;
	double *W;
	int ldw, info;

	if ( n == 1){
		info = LAPACKE_dlarfg_work( m, (&A[0]), &(A[1]), 1, &(T[0]));
	}
	else {

	n1 = n/2;
	n2 = n-n1;

	dgeqr3( m, n1, &(A[0]), lda, T, ldt );
	//dgeqr3_L1( m, n1, &(A[0]), lda, T, ldt );

	/*                                                                                       */
	/* For W we can use the space of T                                                       */
	/* if there is a T computation otherwise we need a buffer                                */
	/*                                                                                       */
	/* Apply Q1^T on A2 with Q = I - Y *T * Y^T                                              */
        /* (see p.613 Table 1 and Table 3 p.615 in Elmroth and Gustavson)                        */
	/* Note that there is a typo in the paper and that (5) is W = Y1 * W.                    */
	/*                                                                                       */
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
	/*                                                                                       */

	ldw=ldt;
	//W = (double *)malloc(ldw*n*sizeof(double)) ; 
	W = &(T[n1*ldt]);

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

	/* from Q1^TA2 to Q2 */
	dgeqr3( m-n1, n2, &(A[n1*lda+n1]), lda, &(T[n1+n1*ldt]), ldt );
	//dgeqr3_L1( m-n1, n2, &(A[n1*m+n1]), lda, &(T[n1+n1*n]), ldt );

	/*                                                                                       */
	/* Computation of the matrix T3 (see p.613 Table 2 in Elmroth and Gustavson)             */
	/*                                                                                       */
	/* ------------------------------------------------------------------------------------- */
	/*     routine  | computation                     | number of floating-point operations  */
	/* ------------------------------------------------------------------------------------- */
	/* (1) DTRMM    | W = Y21^T * Y12                 | n1 * n2^2                            */
	/* (2) DGEMM    | W = W + Y31^T * Y22             | 2 * n1 * n2 * (m-n)                  */
	/* (3) DTRMM    | W = - T1 * W                    | n1^2 * n2                            */
	/* (4) DTRMM    | W = W * T2                      | n1 * n2^2                            */
	/* ------------------------------------------------------------------------------------- */
	/*     total    | W = - T1 * ( Y1^T * Y2 ) * T2   | n1 * n2 * ( 2 * m - n1 )             */
	/* ------------------------------------------------------------------------------------- */
	/*                                                                                       */

	/*      (1) DTRMM    | W = Y21^T * Y12                 | n1 * n2^2                       */

	for (i=0;i<n1;i++) for (j=0;j<n2;j++) T[i+(j+n1)*ldt] = A[j+n1+i*lda];
	cblas_dtrmm ( CblasColMajor, CblasRight, CblasLower, CblasNoTrans, CblasUnit, n1, n2, 1.0e+00, &(A[n1*lda+n1]), lda, &(T[n1*ldt]), ldt); 

	/*      (2) DGEMM    | W = W + Y31^T * Y22             | 2 * n1 * n2 * (m-n)             */

	cblas_dgemm( CblasColMajor, CblasTrans, CblasNoTrans, n1, n2, m-n, +1.0e+00, &(A[n]), lda, &(A[n+n1*lda]),lda, +1.0e+00, &(T[n1*ldt]), ldt);

	/*      (3) DTRMM    | W = - T1 * W                    | n1^2 * n2                       */

	cblas_dtrmm ( CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, n1, n2, -1.0e+00, &(T[0]), ldt, &(T[n1*ldt]), ldt); 

	/*      (4) DTRMM    | W = W * T2                      | n1 * n2^2                       */

	cblas_dtrmm ( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, n1, n2, +1.0e+00, &(T[n1*ldt+n1]), ldt, &(T[n1*ldt]), ldt); 

	}

	return 0;

}
