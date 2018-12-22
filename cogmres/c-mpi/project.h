#include <math.h>
#include <stdlib.h> 
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include "mpi.h"
#include "cblas.h"
#include "lapackwrap.h"
#include "scalapackwrap.h"
#include "PBblacs.h"

extern int cgs_v0                  (int mloc, int n, double *A, int lda, double *R, int ldr, MPI_Comm mpi_comm);
extern int cgs_v1                  (int mloc, int n, double *A, int lda, double *R, int ldr, MPI_Comm mpi_comm);
extern int mgs_v0                  (int mloc, int n, double *A, int lda, double *R, int ldr, MPI_Comm mpi_comm);
extern int mgs_v1                  (int mloc, int n, double *A, int lda, double *R, int ldr, MPI_Comm mpi_comm);
extern int rowmgs_v0               (int mloc, int n, double *A, int lda, double *R, int ldr, MPI_Comm mpi_comm);
extern int rowmgs_v1               (int mloc, int n, double *A, int lda, double *R, int ldr, MPI_Comm mpi_comm);
extern int choleskyqr_A_v0         (int mloc, int n, double *A, int lda, double *R, int ldr, MPI_Comm mpi_comm);
extern int choleskyqr_A_v1         (int mloc, int n, double *A, int lda, double *R, int ldr, MPI_Comm mpi_comm);
extern int choleskyqr_A_v2         (int mloc, int n, double *A, int lda, double *R, int ldr, MPI_Comm mpi_comm);
extern int choleskyqr_A_v3         (int mloc, int n, double *A, int lda, double *R, int ldr, MPI_Comm mpi_comm);
extern int choleskyqr_B_v0         (int mloc, int n, double *A, int lda, double *R, int ldr, MPI_Comm mpi_comm);
extern int choleskyqr_B_v1         (int mloc, int n, double *A, int lda, double *R, int ldr, MPI_Comm mpi_comm);
extern int choleskyqr_B_v2         (int mloc, int n, double *A, int lda, double *R, int ldr, MPI_Comm mpi_comm);
extern int choleskyqr_B_v3         (int mloc, int n, double *A, int lda, double *R, int ldr, MPI_Comm mpi_comm);
extern int reducehouseholder_A_v0  (int mloc, int n, double *A, int lda, double *R, int ldr, MPI_Comm mpi_comm);
extern int reducehouseholder_A_v1  (int mloc, int n, double *A, int lda, double *R, int ldr, MPI_Comm mpi_comm);
extern int reducehouseholder_A_v2  (int mloc, int n, double *A, int lda, double *R, int ldr, MPI_Comm mpi_comm);
extern int reducehouseholder_B_v0  (int mloc, int n, double *A, int lda, double *R, int ldr, MPI_Comm mpi_comm);
extern int reducehouseholder_B_v1  (int mloc, int n, double *A, int lda, double *R, int ldr, MPI_Comm mpi_comm);
extern int reducehouseholder_B_v2  (int mloc, int n, double *A, int lda, double *R, int ldr, MPI_Comm mpi_comm);
extern int reducehouseholder_B_v3  (int mloc, int n, double *A, int lda, double *R, int ldr, MPI_Comm mpi_comm);
extern int scalapackqr2_A          (int mloc, int n, double *A, int lda, double *R, int ldr, MPI_Comm mpi_comm);
extern int scalapackqrf_A          (int mloc, int n, double *A, int lda, double *R, int ldr, MPI_Comm mpi_comm);
extern int scalapackqr2_B          (int mloc, int n, double *A, int lda, double *R, int ldr, MPI_Comm mpi_comm);
extern int scalapackqrf_B          (int mloc, int n, double *A, int lda, double *R, int ldr, MPI_Comm mpi_comm);

extern int check_Rfactor( double *orthlevel, double *represlevel, int mloc, int n, double *A, int lda, double* R, int ldr, MPI_Comm mpi_comm);
extern double check_orthogonality(int mloc, int n, double *Q, int lda, MPI_Comm mpi_comm);
extern double check_orthogonality_RFP(int mloc, int n, double *Q, int lda, MPI_Comm mpi_comm);
extern double check_representativity(int mloc, int n, double *A, int lda, double *Q, int ldq, double* R, int ldr, MPI_Comm mpi_comm);

/* auxiliary routines */
extern int myqrfactorizationofTypeIIImatrices_A ( int n, int p, double *A, double *R);
extern int myqrfactorizationofTypeIIImatrices_B ( int n, int p, double *A, double *R);
extern int LILA_qr_uppers ( int n, double *a, double *b, double *tau, double *work);
extern int LILA_qr_dormqr ( int n, double *tau, double *C2, double *C1, double *work);

extern int LILA_up( int level, int iroot, int mybuddy, int n, double *R, double **H, double **tau, double *work,
	 MPI_Datatype MPI_UPPER, MPI_Comm mpi_comm);
extern int LILA_down( int level, int iroot, int mybuddy, int n, double *R, double *QW, double **H, double **tau, double *work,
	 MPI_Datatype MPI_UPPER, MPI_Comm mpi_comm);

/* MPI Operation for reduction */
extern void LILA_mpiop_sum_upper   ( void *a, void *b, int *len, MPI_Datatype *mpidatatype );
extern void LILA_mpiop_qr_upper_v0 ( void *a, void *b, int *len, MPI_Datatype *mpidatatype );
extern void LILA_mpiop_qr_upper_v1 ( void *a, void *b, int *len, MPI_Datatype *mpidatatype );
extern void LILA_mpiop_qr_upper_v2 ( void *a, void *b, int *len, MPI_Datatype *mpidatatype );
extern void LILA_mpiop_dlapy2      ( void *a, void *b, int *len, MPI_Datatype *mpidatatype );

extern int dgeqr3( int m, int n, double *A, int lda, double *T, int ldt );

