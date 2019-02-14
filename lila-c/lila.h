#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include "cblas.h"
#include "lapacke.h"

extern int dgeqr3( int m, int n, double *A, int lda, double *T, int ldt );

extern int lila_dgeqrf_w02a( int m, int n, int i, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork, int nb );
extern int lila_dgeqrf_w03a( int m, int n, int nb, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );

extern int lila_dgeqrf_w02b( int m, int n, int i, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork, int nb );

extern int lila_dormqrf_w00( int m, int n, int k, int i, int j, int mt, double *A, int lda, double *T, int ldt, double *work, int lwork );
extern int lila_dormqrf_w02( int m, int n, int k, int i, int j, int mt, double *A, int lda, double *T, int ldt, double *work, int lwork );
extern int lila_dormqrf_z00( int m, int n, int k, int i, int j, int mt, double *A, int lda, double *B, int ldb, double *T, int ldt, double *work, int lwork );
extern int lila_dormqrf_z02( int m, int n, int k, int i, int j, int mt, double *A, int lda, double *T, int ldt, double *B, int ldb, double *work, int lwork );

extern int lila_dormqrf_w03( int m, int n, int k, int i, int j, int mt, double *A, int lda, double *T, int ldt, double *work, int lwork );
extern int lila_dormqrf_z03( int m, int n, int k, int i, int j, int mt, double *A, int lda, double *T, int ldt, double *B, int ldb, double *work, int lwork );

extern int lila_wsq_dgeqr2_w02a ( int m, int n,        int i,        int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work );
extern int lila_wsq_dormqrf_z02 ( int m, int n, int k, int i, int j, int mt, double *A, int lda, double *T, int ldt, double *B, int ldb, double *work );
extern int lila_wsq_dormqrbz_w02( int m, int n, int k, int i, int j, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work );

extern int lila_dlarft_connect_w02( int m, int n, int i, int j, int mt, double *A, int lda, double *T, int ldt );
extern int lila_dlarft_connect_w03( int m, int n, int i, int j, int mt, double *A, int lda, double *T, int ldt );

extern int lila_dlarft_w03( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *tau );

extern int lila_dormqrbz_w02( int m, int n, int k, int i, int j, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_dormqrbz_w03( int m, int n, int k, int i, int j, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );

extern int lila_ormhr_w0b( int m, int n, int i, int j, double *A, int lda, double *T, int ldt, double *Q, int ldq, double* *S );
extern int lila_dorgh2( int m, int n, int i, int j, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork, double* S );
extern int lila_ormhr2_3( int m, int n, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_dorgh2_3( int m, int n, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork, double* S );


extern int lila_dgeqrf_w03_level1( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork, int nb );
extern int lila_dorghr_w03_mt( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_dgeqrf_w03_mt_l( int panel, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_dgeqrf_w03_mt( int panel, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_dgeqrf_w03_mt_hr( int panel, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_dgeqr2_w03_l( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_dgeqr2_w03_3( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_dgeqr2_w03_hr( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );

extern int lila_ormhr2_w03_hr( int m, int n, int i, int j, int l, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork, double *S );

extern int lila_dgeqrf_levelx_w03( int panel, int leaf, int n_lvl, int i_lvl, int *nb_lvl, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_dgeqrf_recursive_w03( int panel, int leaf, int nx, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );

extern int lila_wsq_dgeqrf_levelx_w03   ( int panel, int leaf, int n_lvl, int i_lvl, int *nb_lvl, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_wsq_dgeqrf_recursive_w03( int panel, int leaf, int nx, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );

extern int lila_wsq_dormqrbz_w03( int m, int n, int k, int i, int j, int mt, double *A, int lda, double *Q, int ldq, double *T, int ldt, double *work, int lwork );
extern int lila_wsq_dgeqrf_w03_mt_hr( int panel, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_wsq_dormqrf_w03( int m, int n, int k, int i, int j, int mt, double *A, int lda, double *T, int ldt, double *work, int lwork );
extern int lila_wsq_dgeqrf_w03_mt( int panel, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_wsq_dgeqrf_w03_mt_l( int panel, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_wsq_dgeqr2_w03_hr( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_wsq_dgeqr2_w03_3( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_wsq_dgeqr2_w03_l( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_wsq_ormhr2_w03_hr( int m, int n, int i, int j, int l, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork, double *S );



