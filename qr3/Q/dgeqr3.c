#include "qr3.h"

int dgeqr3( int m, int n, double *A, int lda, double *T, int ldt, double *R, int ldr ){

	/* Only the upper part of T is referenced */

	int n1, n2, i, j, info;
	double *A11, *A12, *A21, *A22;
	double *T11, *T12, *T22;
	double *R11, *R12, *R22;

	if ( n == 1){
		info = LAPACKE_dlarfg_work( m, (&A[0]), &(A[1]), 1, &(T[0]));
		(*R) = (*A);
	} else {

	n1 = n/2;
	n2 = n-n1;

	A11 = A;
	A12 = A+n1*lda;
	A21 = A+n1;
	A22 = A+n1*(1+lda);

	T11 = T;
	T12 = T+n1*ldt;
	T22 = T+n1*(1+ldt);

	R11 = R;
	R12 = R+n1*ldt;
	R22 = R+n1*(1+ldt);

	dgeqr3( m, n1, A11, lda, T11, ldt, R11, ldr );


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

	for (i=0;i<n1;i++) for (j=0;j<n2;j++) T12[i+j*ldt] = A12[j*lda+i];
	cblas_dtrmm ( CblasColMajor, CblasLeft, CblasLower, CblasTrans, CblasUnit, n1, n2, 1.0e+00, A11, lda, T12, ldt); 
	/*      (2) DGEMM    | W = W + Y2^T * C2               | 2 * k * ( m - k ) * n           */
	cblas_dgemm( CblasColMajor, CblasTrans, CblasNoTrans, n1, n2, m-n1, +1.0e+00, A21, lda, A22, lda, +1.0e+00, T12, ldt);
	/*      (3) DTRMM    | W = T^T * W                     | k^2 * n                         */
	cblas_dtrmm ( CblasColMajor, CblasLeft, CblasUpper, CblasTrans, CblasNonUnit, n1, n2, 1.0e+00, T11, ldt, T12, ldt); 
	/*      (4) DGEMM    | C2 = C2 - Y2 * W                | 2 * k * ( m - k ) * n           */
	cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, m-n1, n2, n1, -1.0e+00, A21, lda, T12, ldt, +1.0e+00, A22, lda);
	/*      (5) DTRMM    | W = Y1 * W                      | k * ( k - 1 ) * n               */
	cblas_dtrmm ( CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, n1, n2, 1.0e+00, A11, lda, T12, ldt); 
	/*      (6) matrix - | C1 = C1 - W                     | k * n                           */
	for (j=0;j<n2;j++) for (i=0;i<n1;i++) A12[j*lda+i] -= T12[j*ldt+i];
	//for (j=0;j<n2;j++) for (i=0;i<n1;i++) R12[j*lda+i] -= T12[j*ldt+i];
	LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', n1, n2, A12, lda, R12, ldr );

	/* from Q1^TA2 to Q2 */
	dgeqr3( m-n1, n2, A22, lda, T22, ldt, R22, ldr );

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

	/*      (1) DTRMM    | W = Y21^T * Y12                 | n1 * n2^2                       */
	for (i=0;i<n1;i++) for (j=0;j<n2;j++) T12[i+j*ldt] = A21[j+i*lda];
	cblas_dtrmm ( CblasColMajor, CblasRight, CblasLower, CblasNoTrans, CblasUnit, n1, n2, 1.0e+00, A22, lda, T12, ldt); 
	/*      (2) DGEMM    | W = W + Y31^T * Y22             | 2 * n1 * n2 * (m-n)             */
	cblas_dgemm( CblasColMajor, CblasTrans, CblasNoTrans, n1, n2, m-n, +1.0e+00, A21+n2, lda, A22+n2,lda, +1.0e+00, T12, ldt);
	/*      (3) DTRMM    | W = - T1 * W                    | n1^2 * n2                       */
	cblas_dtrmm ( CblasColMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, n1, n2, -1.0e+00, T11, ldt, T12, ldt); 
	/*      (4) DTRMM    | W = W * T2                      | n1 * n2^2                       */
	cblas_dtrmm ( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, n1, n2, +1.0e+00, T22, ldt, T12, ldt); 

	}

	return 0;

}
