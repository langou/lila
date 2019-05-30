#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include "cblas.h"
#include "lapacke.h"

extern int dgeqr3_Q( int m, int n, double *A, int lda, double *T, int ldt );
extern int dgeqr3_Q_ISW( int m, int n, double *A, int lda, double *T, int ldt );




