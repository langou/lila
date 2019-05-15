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

#if !defined(USE_MKL)
extern int dorg2r_( int *m, int *n, int *k, double *a, int *lda, double *tau, double *work, int *info );
//extern void dlarf_( char *side, int *m, int *n, double *v, int *incv, double *tau, double *c, int *ldc, double *work);
#endif

extern int ULTinU( int n, double *L, int ldl, double *U, int ldu );
extern int ApUBTinA( int m, int n, double *A, int lda, double *U, int ldu, double *B, int ldb );
extern int mLUinA( int n, double *A, int lda );

extern int qr3_test_qr_repres_1( double *repres, int m, int n, double *A, int lda, double *Q, int ldq, double *R, int ldr );
extern int qr3_test_qq_orth_1( double *norm_orth_1, int m, int n, double *Q, int ldq  );

extern int our_dorgqr( int m, int n, int k, int nb, double *A, int lda, double *tau, double *work, int lwork );
extern int our_dlarfb_lnfc( int m, int n, int k, double *V, int ldv, double *T, int ldt, double *C, int ldc, double *W );
extern int our_dgeqrf( int m, int n, int nb, double *A, int lda, double *tau, double *work, int lwork );

extern int dV2N( int n, double *T, int ldt );
extern int dN2T( int n, double *tau, double *T, int ldt );
extern int dV2Q( int m, int n, int k, int nb, int nbmin, int nx, double *A, int lda, double *tau, double *work, int lwork, int info );
extern int dVT2Q( int m, int n, double *Q, int ldq  );

extern int dorgqr_after( int m, int n, int k, double *A, int lda, double *T, int ldt, double *Q, int ldq );

extern int dgeqr3R( int m, int n, double *A, int lda, double *T, int ldt, double *R, int ldr );
extern int dgeqr3_ISW( int m, int n, double *A, int lda, double *T, int ldt );
extern int dgeqr3( int m, int n, double *A, int lda, double *T, int ldt );

extern int qr3_larft( int m, int n, double *A, int lda, double *T, int ldt, double *tau ); 
extern int qr3_dorgqr( int m, int n, double *A, int lda, double *T, int ldt, double *tau );

extern int qr3_dA2QRTV_fake( int m, int n, double *A, int lda, double *Q, int ldq, double *T, int ldt );
extern int qr3_dA2QRTV( int m, int n, double *A, int lda, double *Q, int ldq, double *T, int ldt );

extern int qr3_null_dorgqr( int m, int n, int k, int nb, double *A, int lda, double *Q, int ldq, double *tau, double *work, int lwork );
extern int qr3_dorgqr_level1( int m, int n, int k, int nb, double *A, int lda, double *tau, double *work, int lwork );
extern int qr3_dorgqr_level1_UT( int m, int n, int k, int nb, double *A, int lda, double *tau, double *work, int lwork );
extern int UinvLTinU( int n, double *L, int ldl, double *U, int ldu );
extern int UinvLTinU_cheat( int n, double *L, int ldl, double *U, int ldu );
extern int dVS2Q( int m, int n, double *Q, int ldq  );

extern int our_dorgqr_Q2( int m, int n, int k, int nb, double *A, int lda, double *tau, double *work, int lwork );
extern int our_dorg2r_Q2( int m, int n, int k, double *A, int lda, double *tau, double *work, int lwork );
extern int our_dorg2r( int m, int n, int k, double *A, int lda, double *tau, double *work, int lwork );

extern int wrapper_dlarf( char side, int m, int n, double *V, int incv, double tau, double *C, int ldc, double *work);
