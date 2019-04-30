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

extern long int flops_org2r( int m, int n, int k );
extern long int flops_gemm( int m, int n, int k );
extern long int flops_trmm( char S, int m, int n );
extern long int flops_syrk( int n, int k );
extern long int flops_larft( int m, int k );
extern long int flops_larfb( int m, int n, int k );

extern long int flops_mLUinA( int n );
extern long int flops_V2N( int n );
extern long int flops_ApUBTinA( int m, int n );
extern long int flops_ULTinU( int n );
extern long int flops_N2T( int n );
extern long int flops_VT2Q( int m, int n );

extern long int flops_mLUinA_check( int n );
extern long int flops_V2N_check( int n );
extern long int flops_ApUBTinA_check( int m, int n );
extern long int flops_ULTinU_check( int n );
extern long int flops_N2T_check( int n );
extern long int flops_VT2Q_check( int m, int n );

extern long int flops_lapack_orgqr( int m, int n, int k, int nb );
extern long int flops_lapack_orgqr_check( int m, int n, int k, int nb );

extern long int flops_lapack_larfb( int m, int n, int k );
extern long int flops_lapack_larfb_check( int m, int n, int k );

extern long int flops_lapack_orgqr_from_last_org2r( int int_m, int int_n, int int_k, int int_b );
extern long int flops_lapack_orgqr_from_first_org2r( int int_m, int int_n, int int_k, int int_b );
extern long int flops_lapack_orgqr_from_larft( int int_m, int int_n, int int_k, int int_b );

extern long int flops_org2r_n1( int int_m );
extern long int flops_lapack_larf( int int_m, int int_n );
extern long int flops_lapack_larfg( int int_m );
extern long int flops_org2r_check( int m, int n, int k );

extern long int flops_lapack_geqrf_check( int m, int n, int nb );
extern long int flops_lapack_geqr2( int int_m, int int_n );
extern long int flops_lapack_geqrf( int int_m, int int_n, int int_nb );

extern long int flops_dorgqr_after_check( int m, int n, int k );
extern long int flops_dorgqr_after( int m, int n, int k );
 
extern int dgeqr3_right( int m, int n, double *A, int lda, double *T, int ldt );
extern int dgeqr3_left( int m, int n, double *A, int lda, double *T, int ldt );

extern int qr3_larft( int m, int n, double *A, int lda, double *T, int ldt, double *tau ); 
extern int qr3_dorgqr( int m, int n, double *A, int lda, double *T, int ldt, double *work, int lwork, double *tau );



