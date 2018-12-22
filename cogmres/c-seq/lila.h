#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include "cblas.h"
#include "lapacke.h"

extern int matvec_A ( int n, double *y, double *x );
extern int matvec_Ml( int n, double *y, double *x );
extern int matvec_Mr( int n, double *y, double *x );
extern int dgmres( int n, double *b, double *x, int m, int max_it, double tol );
