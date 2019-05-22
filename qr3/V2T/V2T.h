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

extern int dV2N( int n, double *T, int ldt );
extern int dN2T( int n, double *tau, double *T, int ldt );

extern int qr3_larft    ( int m, int n, double *A, int lda, double *T, int ldt, double *tau );
extern int qr3_larft_ISW( int m, int n, double *A, int lda, double *T, int ldt, double *tau );

extern int qr3_larft_UT    ( int m, int n, double *A, int lda, double *T, int ldt, double *tau );
extern int qr3_larft_UT_ISW( int m, int n, double *A, int lda, double *T, int ldt, double *tau );

extern int qr3_dorgqr   ( int m, int n, double *Q, int ldq, double *T, int ldt, double *tau );
extern int qr3_dorgqr_UT( int m, int n, double *Q, int ldq, double *T, int ldt, double *tau );

extern int V2T_test_qq_orth_1  ( double *norm_orth_1, int m, int n, double *Q, int ldq  );
extern int V2T_test_qr_repres_1( double *norm_repres, int m, int n, double *A, int lda, double *Q, int ldq, double *R, int ldr );

extern int dVS2Q    ( int m, int n, double *Q, int ldq  );
extern int mLUinA   ( int n, double *A, int lda );
extern int UinvLTinU( int n, double *L, int ldl, double *U, int ldu );

extern int qr3_larft_ISW_V2T   ( int m, int n, double *A, int lda, double *T, int ldt, double *tau );
extern int qr3_larft_ISW_V2T_UT( int m, int n, double *A, int lda, double *T, int ldt, double *tau );



