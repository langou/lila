#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include "cblas.h"
#include "lapacke.h"

extern int dgeqr3( int m, int n, double *A, int lda, double *T, int ldt );

extern int lila_dgeqr2_w02a( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );

extern int lila_dormqrf_w00( int m, int n, int k, int i, int j, int mt, double *A, int lda, double *T, int ldt, double *work, int lwork );
extern int lila_dormqrf_w02( int m, int n, int k, int i, int j, int mt, double *A, int lda, double *T, int ldt, double *work, int lwork );
extern int lila_dormqrf_z00( int m, int n, int k, int i, int j, int mt, double *A, int lda, double *B, int ldb, double *T, int ldt, double *work, int lwork );
extern int lila_dormqrf_z02( int m, int n, int k, int i, int j, int mt, double *A, int lda, double *T, int ldt, double *B, int ldb, double *work, int lwork );

extern int  lila_wsq_dgeqr2_w02a( int m, int n,        int i,        int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work );
extern int  lila_wsq_dormqrf_z02( int m, int n, int k, int i, int j, int mt, double *A, int lda, double *T, int ldt, double *B, int ldb, double *work );
extern int lila_wsq_dormqrbz_w02( int m, int n, int k, int i, int j, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work );


extern int lila_dlarft_connect_w02( int m, int n, int i, int j, int mt, double *A, int lda, double *T, int ldt );

extern int lila_dormqrbz_w02( int m, int n, int k, int i, int j, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );



