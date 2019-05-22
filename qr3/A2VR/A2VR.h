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

extern int dgeqr3R_ISW( int m, int n, double *A, int lda, double *T, int ldt, double *R, int ldr );
extern int dgeqr3R    ( int m, int n, double *A, int lda, double *T, int ldt, double *R, int ldr );
extern int dgeqr3R_UT_ISW( int m, int n, double *A, int lda, double *T, int ldt, double *R, int ldr );
extern int dgeqr3R_UT    ( int m, int n, double *A, int lda, double *T, int ldt, double *R, int ldr );

extern int qr3_dorgqr( int m, int n, double *Q, int ldq, double *T, int ldt, double *tau );

extern int A2VR_test_qq_orth_1( double *norm_orth_1, int m, int n, double *Q, int ldq  );
extern int A2VR_test_qr_repres_1( double *norm_repres, int m, int n, double *A, int lda, double *Q, int ldq, double *R, int ldr );

extern int dV2tau( int m, int n, double *A, int lda, double *tau );



