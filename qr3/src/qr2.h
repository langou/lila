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

extern int qr2_dgeqr3R( int m, int n, double *A, int lda, double *T, int ldt, double *R, int ldr );
extern int qr2_dgeqr3R_ISW( int m, int n, double *A, int lda, double *T, int ldt, double *R, int ldr );
