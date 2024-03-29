#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
//#include "cblas.h"
//#include "lapacke.h"
#include "mkl_cblas.h"
#include "mkl_lapacke.h"

extern int dgeqr3( int m, int n, double *A, int lda, double *T, int ldt );

extern int lila_dgetrf_b( int m, int n, int i, int mt, double *A, int lda, double *Q, int ldq, double *work, int lwork );

extern int lila_dormqrf_w00( int m, int n, int k, int i, int j, int mt, double *A, int lda, double *T, int ldt, double *work, int lwork );
extern int lila_dormqrf_w02( int m, int n, int k, int i, int j, int mt, double *A, int lda, double *T, int ldt, double *work, int lwork );
extern int lila_dormqrf_w03( int m, int n, int k, int i, int j, int mt, double *A, int lda, double *T, int ldt, double *work, int lwork );

extern int lila_dormqrf_z00( int m, int n, int k, int i, int j, int mt, double *A, int lda, double *B, int ldb, double *T, int ldt, double *work, int lwork );
extern int lila_dormqrf_z02( int m, int n, int k, int i, int j, int mt, double *A, int lda, double *B, int ldb, double *T, int ldt, double *work, int lwork );
extern int lila_dormqrf_z03( int m, int n, int k, int i, int j, int mt, double *A, int lda, double *B, int ldb, double *T, int ldt, double *work, int lwork );

extern int lila_dgeqrf_w03_mt( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_ormhr2_w03_hr( int m, int n, int i, int j, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork, int *S );

extern int lila_dgeqr1_w03a( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_dgeqr1_w03b( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );

extern int lila_dgeqr2_w02a( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_dgeqr2_w02b( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );

extern int lila_dgeqr2_w03a( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_dgeqr2_w03b( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );

extern int lila_dlarft_connect_w02( int m, int n, int i, int j, int mt, double *A, int lda, double *T, int ldt );
extern int lila_dlarft_connect_w03( int m, int n, int i, int j, int mt, double *A, int lda, double *T, int ldt );

extern int lila_dlarft_w03( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *tau );
extern int lila_dlarft_w03_b( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt );

extern int lila_ormhr_w0b( int m, int n, int i, int j, double *A, int lda, double *T, int ldt, double *Q, int ldq, int *S );
extern int lila_dorgh2( int m, int n, int j, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork, int *S );
extern int lila_dorghr( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork, int *S );

extern int lila_dormqrbz_w00( int m, int n, int k, int i, int j, int mt, double *A, int lda, double *Q, int ldq, double *T, int ldt, double *work, int lwork );
extern int lila_dormqrbz_w02( int m, int n, int k, int i, int j, int mt, double *A, int lda, double *Q, int ldq, double *T, int ldt, double *work, int lwork );
extern int lila_dormqrbz_w03( int m, int n, int k, int i, int j, int mt, double *A, int lda, double *Q, int ldq, double *T, int ldt, double *work, int lwork );

extern int lila_dgeqrf_levelx_w00( int n_lvl, int i_lvl, int *nb_lvl, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_dgeqrf_levelx_w02( int n_lvl, int i_lvl, int *nb_lvl, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_dgeqrf_levelx_w03( int n_lvl, int i_lvl, int *nb_lvl, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );

extern int lila_dgeqrf_recursive_w02( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_dgeqrf_recursive_w03( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
