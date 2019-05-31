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

extern int check_qq_orth( double *norm_orth_1, int m, int n, double *Q, int ldq  );
extern int check_qr_repres( double *norm_repres, int m, int n, double *A, int lda, double *Q, int ldq, double *R, int ldr );
extern int check_q2A_repres( double *norm_repres, int m, int n, int k, double *A, int lda, double *Q2, int ldq2 );




