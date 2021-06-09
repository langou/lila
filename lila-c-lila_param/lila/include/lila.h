#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>


#include "cblas.h"
#include "lapacke.h"
//#include "mkl_cblas.h"
//#include "mkl_lapacke.h"

extern int dgeqr3( int m, int n, double *A, int lda, double *T, int ldt );

extern int lila_dormqrf_z00( int m, int n, int k, int i, int j, int mt, double *A, int lda, double *B, int ldb, double *T, int ldt, double *work, int lwork );
extern int lila_dormqrf_z02( int m, int n, int k, int i, int j, int mt, double *A, int lda, double *T, int ldt, double *B, int ldb, double *work, int lwork );
extern int lila_dormqrf_z03( int m, int n, int k, int i, int j, int mt, double *A, int lda, double *T, int ldt, double *B, int ldb, double *work, int lwork );
extern int lila_dormqrf_w00( int m, int n, int k, int i, int j, int mt, double *A, int lda, double *T, int ldt, double *work, int lwork );
extern int lila_dormqrf_w03( int *lila_param, int m, int n, int k, int i, int j, int mt, double *A, int lda, double *T, int ldt, double *work, int lwork );

extern int lila_dormqrbz_w03( int *lila_param, int m, int n, int k, int i, int j, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );

extern int lila_ormhr2_3 ( int m, int n, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_ormhr_w0b( int m, int n, int i, int j, double *A, int lda, double *T, int ldt, double *Q, int ldq, double* *S );
extern int lila_dorgh2   ( int m, int n, int i, int j, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork, double* S );
extern int lila_dorgh2_3 ( int m, int n, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork, double* S );
extern int lila_ormhr2_w03_hr( int m, int n, int i, int j, int l, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork, double *S );


extern int lila_query_dgeqrf_w03_levelx   ( int *lila_param, int n_lvl, int i_lvl, int *nb_lvl, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_query_dgeqrf_w03_recursive( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_query_dgeqrf_ker_recursive( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_query_dgeqrf_w03_mt       ( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_query_dgeqrf_w03_mt_hh    ( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_query_dormqrbz_w03        ( int m, int n, int k, int i, int j, int mt, double *A, int lda, double *Q, int ldq, double *T, int ldt, double *work, int lwork ); 
extern int lila_query_dgeqrf_w03_mt_hr    ( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_query_dgeqrf_w03_mt_l     ( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_query_dgeqr2_w03_hr       ( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_query_dgeqr2_w03_3        ( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_query_dgeqr2_w03_l        ( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_query_ormhr2_w03_hr       ( int m, int n, int i, int j, int l, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork, double *S );
extern int lila_query_dorgh2_3            ( int m, int n, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork, double *S );
extern int lila_query_dormqrf_w03         ( int m, int n, int k, int i, int j, int mt, double *A, int lda, double *T, int ldt, double *work, int lwork );
extern int lila_query_dgeqrf_LAPACK_appendcols( int m, int k, int n, int mt, double *A, int lda, double *Q, int ldq, double *T, int ldt, double *work, int lwork );
/*

extern int lila_dlarft_connect_w03( int *lila_param, int m, int n, int i, int j, int mt, double *A, int lda, double *T, int ldt );
extern int lila_dlarft_w03        ( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *tau ); 

extern int lila_dgeqrf_w03_levelx   ( int *lila_param, int n_lvl, int i_lvl, int *nb_lvl, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_dgeqrf_w03_recursive( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_dgeqr2_w03          ( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_dgeqrf_w03_mt       ( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_dgeqrf_w03_mt_hh    ( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_dorghr_w03_mt                        ( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_dgeqrf_w03_mt_l     ( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_dgeqrf_w03_mt_hr    ( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_dgeqr2_w03_l        ( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_dgeqr2_w03_3        ( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_dgeqr2_w03_hr       ( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );

extern int lila_dgeqrf_v03_levelx   ( int *lila_param, int n_lvl, int i_lvl, int *nb_lvl, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *work, int lwork );
extern int lila_dgeqrf_v03_recursive( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *work, int lwork );
extern int lila_dgeqrf_v03_mt       ( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *work, int lwork );
extern int lila_dgeqr2_v03          ( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *work, int lwork );
extern int lila_dgeqrf_v03_mt_hh    ( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *work, int lwork );
extern int lila_dgeqrf_v03_mt_l     ( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *work, int lwork );
extern int lila_dgeqr2_v03_l        ( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *work, int lwork );
extern int lila_dgeqr2_v03_3        ( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *work, int lwork );

extern int lila_dgeqrf_q03_levelx   ( int *lila_param, int n_lvl, int i_lvl, int *nb_lvl, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_dgeqrf_q03_recursive( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_dgeqrf_q03_mt_hh    ( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_dgeqrf_q03_mt       ( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_dgeqrf_q03_mt_l     ( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_dgeqr2_q03          ( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_dgeqr2_q03_l        ( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_dgeqr2_q03_3        ( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );

extern int lila_dgeqrf_t03  ( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *work, int lwork );
extern int lila_dgeqr2_t03_l( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *tau, double *work, int lwork );
extern int lila_dgeqrf_t03_l( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *tau, double *work, int lwork );
*/

extern int lila_dT2tau_w03( int m, int n, int i, int mt, double *T, int ldt, double *tau );
extern int lila_dV2tau_w03( int m, int n, int i, int mt, double *A, int lda, double *tau );

extern int lila_dgeqrf_ker_recursive( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );

extern int lila_dgeqrf_qr2_recursive2( int *lila_param, int m, int n, int i, int mt, int nb1, int nb2, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_dgeqrf_qr2_recursive ( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_dgeqrf_qr2           ( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );

extern int lila_dgeqrf_w03_appendcols   ( int *lila_param, int m, int k, int n, int mt, double *A, int lda, double *Q, int ldq, double *T, int ldt, double *work, int lwork );
extern int lila_dgeqrf_LAPACK_appendcols( int m, int k, int n, int mt, double *A, int lda, double *Q, int ldq, double *T, int ldt, double *work, int lwork );

extern int lila_dgeqrf_appendcols       ( int *lila_param, int m, int n, int k, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );


extern int lila_dgeqrf_recursive( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_dgeqrf_levelx   ( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_dgeqrf          ( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_dgeqrf_mt       ( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_dgeqrf_mt_l     ( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_dgeqrf_mt_hr    ( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_dgeqrf_mt_hh    ( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_dgeqr2          ( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_dgeqr2_l        ( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_dgeqr2_3        ( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_dgeqr2_hr       ( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );
extern int lila_dlarft_connect  ( int *lila_param, int m, int n, int i, int j, int mt, double *A, int lda, double *T, int ldt );
extern int lila_dormqrbz        ( int *lila_param, int m, int n, int k, int i, int j, int mt, double *A, int lda, double *Q, int ldq, double *T, int ldt, double *work, int lwork );
extern int lila_dormqrf         ( int *lila_param, int m, int n, int k, int i, int j, int mt, double *A, int lda, double *T, int ldt, double *work, int lwork );
extern int lila_dlarft          ( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *tau );
extern int lila_dT2tau          ( int *lila_param, int m, int n, int i, int mt, double *T, int ldt, double *tau );

extern double lila_test_qr_repres_1( int m, int n, int i, double *A, int lda, double *Q, int ldq, double *R, int ldr );
extern double lila_test_qq_orth_1  ( int m, int n, int i, double *Q, int ldq );
extern double lila_test_r_repres_2 ( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *R, int ldr );
extern double lila_test_hh_repres  ( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *V, int ldv );
extern double lila_test_vt_repres  ( int *lila_param, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *V, int ldv );

