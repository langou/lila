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

extern int qr2_dgeqr3R       ( int m, int n, double *A, int lda, double *T, int ldt, double *R, int ldr );
extern int qr2_dgeqr3R_UT    ( int m, int n, double *A, int lda, double *T, int ldt, double *R, int ldr );
extern int qr2_dgeqr3R_ISW   ( int m, int n, double *A, int lda, double *T, int ldt, double *R, int ldr );
extern int qr2_dgeqr3R_UT_ISW( int m, int n, double *A, int lda, double *T, int ldt, double *R, int ldr );

extern int qr2_dorgqr   ( int m, int n, double *Q, int ldq, double *T, int ldt, double *tau );
extern int qr2_dorgqr_UT( int m, int n, double *Q, int ldq, double *T, int ldt, double *tau );

extern int qr2_dV2tau( int m, int n, double *A, int lda, double *tau );

extern int qr2_larft           ( int m, int n, double *A, int lda, double *T, int ldt, double *tau );
extern int qr2_larft_UT        ( int m, int n, double *A, int lda, double *T, int ldt, double *tau );
extern int qr2_larft_ISW       ( int m, int n, double *A, int lda, double *T, int ldt, double *tau );
extern int qr2_larft_UT_ISW    ( int m, int n, double *A, int lda, double *T, int ldt, double *tau );
extern int qr2_larft_ISW_V2T   ( int m, int n, double *A, int lda, double *T, int ldt, double *tau );
extern int qr2_larft_ISW_V2T_UT( int m, int n, double *A, int lda, double *T, int ldt, double *tau );

extern int qr2_dV2N     ( int n, double *T, int ldt );
extern int qr2_dN2T     ( int n, double *tau, double *T, int ldt );
extern int qr2_dVS2Q    ( int m, int n, double *Q, int ldq  );
extern int qr2_mLUinA   ( int n, double *A, int lda );
extern int qr2_UinvLTinU( int n, double *L, int ldl, double *U, int ldu );







