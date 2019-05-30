#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include "cblas.h"
#include "lapacke.h"



extern int qr3_aux_dorgq2r( int m, int n, int k, double *A, int lda, double *T, int ldt, double *Q, int ldq );
extern int lapack_our_dorgq2r( int m, int n, int k, int nb, double *A, int lda, double *Q, int ldq, double *tau, double *work, int lwork );




