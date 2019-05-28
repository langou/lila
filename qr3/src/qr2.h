#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>

#if !defined(USE_MKL)
#include "cblas.h"
#include "lapacke.h"
extern void dlarf_( char *side, int *m, int *n, double *v, int *incv, double *tau, double *c, int *ldc, double *work);
extern void dorg2r_( int *m, int *n, int *k, double *a, int *lda, double *tau, double *work, int *info );
#endif

#if defined(USE_MKL)
#include "mkl.h"
#endif

extern int qr2_dgeqr3_R       ( int m, int n, double *A, int lda, double *T, int ldt, double *R, int ldr );
extern int qr2_dgeqr3_R_UT    ( int m, int n, double *A, int lda, double *T, int ldt, double *R, int ldr );
extern int qr2_dgeqr3_R_ISW   ( int m, int n, double *A, int lda, double *T, int ldt, double *R, int ldr );
extern int qr2_dgeqr3_R_UT_ISW( int m, int n, double *A, int lda, double *T, int ldt, double *R, int ldr );

extern int qr2_dorgqr3     ( int m, int n, double *Q, int ldq, double *T, int ldt, double *tau );
extern int qr2_dorgqr3_UT  ( int m, int n, double *Q, int ldq, double *T, int ldt, double *tau );
extern int qr2_dorgqr3_VT2Q( int m, int n, double *Q, int ldq, double *T, int ldt, double *tau );
extern int qr2_dorgqr3_VS2Q( int m, int n, double *Q, int ldq, double *T, int ldt, double *tau );

extern int qr2_dV2tau( int m, int n, double *A, int lda, double *tau );

extern int qr2_larft3           ( int m, int n, double *A, int lda, double *T, int ldt, double *tau );
extern int qr2_larft3_UT        ( int m, int n, double *A, int lda, double *T, int ldt, double *tau );
extern int qr2_larft3_ISW       ( int m, int n, double *A, int lda, double *T, int ldt, double *tau );
extern int qr2_larft3_UT_ISW    ( int m, int n, double *A, int lda, double *T, int ldt, double *tau );
extern int qr2_larft3_ISW_V2T   ( int m, int n, double *A, int lda, double *T, int ldt, double *tau );
extern int qr2_larft3_ISW_V2T_UT( int m, int n, double *A, int lda, double *T, int ldt, double *tau );

extern int qr2_dV2N       ( int n, double *T, int ldt );
extern int qr2_dN2T       ( int n, double *tau, double *T, int ldt );
extern int qr2_dVS2Q      ( int m, int n, double *Q, int ldq  );
extern int qr2_dVT2Q      ( int m, int n, double *Q, int ldq  );
extern int qr2_mLUinA     ( int n, double *A, int lda );
extern int qr2_ULTinU     ( int n, double *L, int ldl, double *U, int ldu );
extern int qr2_ApUBTinA   ( int m, int n, double *A, int lda, double *U, int ldu, double *B, int ldb );
extern int qr2_UinvLTinU  ( int n, double *L, int ldl, double *U, int ldu );

extern int level1_lapack_dorgqr_Q2( int m, int n, int k, int nb, double *A, int lda, double *tau, double *work, int lwork );

extern int our_dlarfb_lnfc( int m, int n, int k, double *V, int ldv, double *T, int ldt, double *C, int ldc, double *W );
extern int our_dorg2r_Q2( int m, int n, int k, double *A, int lda, double *tau, double *work, int lwork );
extern int wrapper_dlarf( char side, int m, int n, double *V, int incv, double tau, double *C, int ldc, double *work);









