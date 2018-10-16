#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include "cblas.h"
#include "lapacke.h"

extern int dgeqr3( int m, int n, double *A, int lda, double *T, int ldt );

extern int lila_dge_qr_wq_vL0( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );

extern int lila_dge_qr_ormqrf_vL0( int m, int n, int k, int i, int j, int mt, double *A, int lda, double *T, int ldt, double *work, int lwork );

extern int lila_dge_qr_ormqrbz_vL0( int m, int n, int k, int i, int j, int mt, double *A, int lda, double *Q, int ldq, double *T, int ldt, double *work, int lwork );

extern int lila_dge_qr_vr0( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq );

extern int lila_dge_qr_wq_manylevels_vr0( int n_lvl, int i_lvl, int *nb_lvl, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq );

extern int lila_dge_qr_wq_manylevels_INTERCEPT_level1_vr0( int n_lvl, int i_lvl, int *nb_lvl, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );

extern int lila_dge_qr_wq_levelx_w00( int n_lvl, int i_lvl, int *nb_lvl, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *Q, int ldq, double *work, int lwork );

extern int lila_dge_qr_ormqrf_w02( int m, int n, int k, int i, int j, int mt, double *A, int lda, double *T, int ldt, double *work, int lwork );

extern int lila_dge_qr_ormqrf_w03( int m, int n, int k, int i, int j, int mt, double *A, int lda, double *T, int ldt, double *TTT, int llldddttt, double *work, int lwork );

extern int lila_dge_qr_ormqrbz_w02( int m, int n, int k, int i, int j, int mt, double *A, int lda, double *Q, int ldq, double *T, int ldt, double *work, int lwork );

extern int lila_dge_qr_wq_w02( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *TTT, int llldddttt, double *Q, int ldq, double *work, int lwork );

extern int lila_dge_qr_larft_connect_w02( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt );

extern int lila_dge_qr_larft_connect_w03( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt );

extern int lila_dge_qr_ormqrf_w00( int m, int n, int k, int i, int j, int mt, double *A, int lda, double *T, int ldt, double *work, int lwork );

extern int lila_dge_qr_ormqrbz_w00( int m, int n, int k, int i, int j, int mt, double *A, int lda, double *Q, int ldq, double *T, int ldt, double *work, int lwork );

extern int lila_dge_qr_wq_levelx_w03( int n_lvl, int i_lvl, int *nb_lvl, int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *TTT, int llldddttt, double *Q, int ldq, double *work, int lwork );

extern int lila_dlarft_w03( int m, int n, int i, int mt, double *A, int lda, double *T, int ldt, double *tau);

