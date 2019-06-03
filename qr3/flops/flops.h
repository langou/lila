#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>

extern long int flops_gemm( int m, int n, int k );
extern long int flops_trmm( char S, int m, int n );
extern long int flops_syrk( int n, int k );
extern long int flops_larft( int m, int k );
extern long int flops_larfb( int m, int n, int k );

extern long int flops_mLUinA( int n );
extern long int flops_V2N( int n );
extern long int flops_ApUBTinA( int m, int n );
extern long int flops_ULTinU( int n );
extern long int flops_N2T( int n );
extern long int flops_VT2Q( int m, int n );

extern long int flops_mLUinA_check( int n );
extern long int flops_V2N_check( int n );
extern long int flops_ApUBTinA_check( int m, int n );
extern long int flops_ULTinU_check( int n );
extern long int flops_N2T_check( int n );
extern long int flops_VT2Q_check( int m, int n );

extern long int flops_lapack_larfb( int m, int n, int k );
extern long int flops_lapack_larfb_check( int m, int n, int k );

extern long int flops_lapack_larf( int m, int n );
extern long int flops_lapack_larfg( int m );

extern long int flops_lapack_geqr2_check( int m, int n );
extern long int flops_lapack_geqr2( int m, int n );
extern long int flops_lapack_geqrf_check( int m, int n, int nb );
extern long int flops_lapack_geqrf( int m, int n, int nb );

extern long int flops_lapack_orgqr_check( int m, int n, int k, int nb );
extern long int flops_lapack_org2r_check( int m, int n, int k );

extern long int flops_lapack_org2r( int int_m, int int_n, int int_k );

extern long int flops_legacy_lapack_org2r( int m, int n, int k );
extern long int flops_legacy_lapack_org2r_check( int m, int n, int k );
extern long int flops_legacy_lapack_larf( int m, int n );
extern long int flops_legacy_lapack_geqr2_check( int m, int n );
extern long int flops_legacy_lapack_geqr2( int m, int n );

extern long int flops_legacy_lapack_larfg( int int_m );

extern long int flops_lapack_orgqr_bef_check( int m, int n, int k, int nb );
extern long int flops_lapack_orgqr_bef( int m, int n, int k, int nb );

extern long int flops_lapack_org2r_Q1_check( int m, int n, int k );
extern long int flops_lapack_org2r_Q1( int m, int n, int k );
extern long int flops_lapack_org2r_Q2_check( int m, int n, int k );
extern long int flops_lapack_org2r_Q2( int m, int n, int k );

extern long int flops_lapack_geqrf_bef( int m, int n, int nb );
extern long int flops_lapack_geqrf_bef_check( int m, int n, int nb );

extern long int flops_geqr3( int m, int n );
extern long int flops_geqr3_check( int m, int n );
extern long int flops_geqr3_wob_check( int m, int n );
extern long int flops_geqr3_bef_constructT_check( int m, int n );
extern long int flops_geqr3_bef_useT_check( int n );
extern long int flops_geqr3_bef_constructT( int m, int n );
extern long int flops_geqr3_bef_useT( int n );
extern long int flops_geqr3_noR( int int_m, int int_n );
extern long int flops_geqr3_noR_check( int int_m, int int_n );
extern long int flops_geqr3_onlyR( int int_n );
extern long int flops_geqr3_onlyR_check( int int_n );

extern long int flops_geqr3_ISW( int m, int n );
extern long int flops_geqr3_ISW_check( int m, int n );
extern long int flops_geqr3_ISW_save( int m, int n );
extern long int flops_geqr3_ISW_constructT_check( int m, int n );
extern long int flops_geqr3_ISW_constructT( int int_m, int int_n );

extern long int flops_geqr3_UT_check( int int_m, int int_n );
extern long int flops_geqr3_UT_save_check( int int_m, int int_n );
extern long int flops_geqr3_UT( int int_m, int int_n );

extern long int flops_geqr3_ISW_UT_check( int int_m, int int_n );
extern long int flops_geqr3_ISW_UT( int int_m, int int_n );

extern long int flops_qr2_dorgqr_check( int m, int n );
extern long int flops_qr2_dorgqr( int int_m, int int_n );

extern long int flops_dlarft3_check( int int_m, int int_n );
extern long int flops_dlarft3( int int_m, int int_n );





