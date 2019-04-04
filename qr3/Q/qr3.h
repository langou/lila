#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include "cblas.h"
#include "lapacke.h"

extern int ULTinU( int n, double *L, int ldl, double *U, int ldu );
extern int ApUBTinA( int m, int n, double *A, int lda, double *U, int ldu, double *B, int ldb );
extern int mLUinA( int n, double *A, int lda );

extern int qr3_test_qr_repres_1( double *repres, int m, int n, double *A, int lda, double *Q, int ldq, double *R, int ldr );
extern int qr3_test_qq_orth_1( double *norm_orth_1, int m, int n, double *Q, int ldq  );

extern int dorgqr( int m, int n, int k, int nb, int nbmin, int nx, double *A, int lda, double *tau, double *work, int lwork, int info );
extern int xV2N( int n, double *T, int ldt );
extern int xN2T( int n, double *tau, double *T, int ldt );



extern unsigned long int flops_orgqr( int m, int n, int k );
extern unsigned long int flops_gemm( int m, int n, int k );
extern unsigned long int flops_trmm( int m, int n, char S );
extern unsigned long int flops_syrk( int n, int k );

