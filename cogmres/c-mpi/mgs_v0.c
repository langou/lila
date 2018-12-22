#include "project.h"

int mgs_v0(int mloc, int n, double *A, int lda, double *R, int ldr, MPI_Comm mpi_comm){

	int i,j;

	for ( j = 0; j < n; j++) {
		for ( i = 0; i < j; i++) {
			R[j*n+i] = cblas_ddot(mloc,&(A[i*mloc]),1,&(A[j*mloc]),1) ;
			MPI_Allreduce( MPI_IN_PLACE, &(R[j*n+i]), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			cblas_daxpy(mloc,-R[j*n+i],&(A[i*mloc]),1,&(A[j*mloc]),1);
		}
		R[j*n+j] = cblas_ddot(mloc,&(A[j*mloc]),1,&(A[j*mloc]),1);
		MPI_Allreduce( MPI_IN_PLACE, &(R[j*n+j]), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		R[j*n+j] = sqrt(R[j*n+j]);
		cblas_dscal(mloc,1/R[j*n+j],&(A[j*mloc]),1);
	}

	return 0;

}
