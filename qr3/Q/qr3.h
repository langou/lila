#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>

#if !defined(USE_MKL)
#include "cblas.h"
#include "lapacke.h"
#endif

#if defined(USE_MKL)
#include "mkl.h"
#endif

//extern int dorg2r_( int *m, int *n, int *k, double *a, int *lda, double *tau, double *work, int *info );

extern int ULTinU( int n, double *L, int ldl, double *U, int ldu );
extern int ApUBTinA( int m, int n, double *A, int lda, double *U, int ldu, double *B, int ldb );
extern int mLUinA( int n, double *A, int lda );

extern int qr3_test_qr_repres_1( double *repres, int m, int n, double *A, int lda, double *Q, int ldq, double *R, int ldr );
extern int qr3_test_qq_orth_1( double *norm_orth_1, int m, int n, double *Q, int ldq  );

extern int our_dorgqr( int m, int n, int k, int nb, int nbmin, int nx, double *A, int lda, double *tau, double *work, int lwork, int info );
extern int dV2N( int n, double *T, int ldt );
extern int dN2T( int n, double *tau, double *T, int ldt );
extern int dV2Q( int m, int n, int k, int nb, int nbmin, int nx, double *A, int lda, double *tau, double *work, int lwork, int info );
extern int dVT2Q( int m, int n, double *Q, int ldq  );

extern unsigned long int flops_org2r( int m, int n, int k );
extern unsigned long int flops_gemm( int m, int n, int k );
extern unsigned long int flops_trmm( int m, int n, char S );
extern unsigned long int flops_syrk( int n, int k );
extern unsigned long int flops_larft( int m, int k );
extern unsigned long int flops_larfb( int m, int n, int k );

extern unsigned long int flops_V2N( int n );
extern unsigned long int flops_VT2Q( int m, int n );
extern unsigned long int flops_N2T( int n );



