#include "project.h"

/* This version will overflow if a_{i,j} > sqrt(HUGE) */

double check_representativity(int mloc, int n, double *A, int lda,
	 double *Q, int ldq, double* R, int ldr, MPI_Comm mpi_comm){

	/* A, Q and R are input only (unchanged on exit) */
	/* Data distribution of the matrix Q is by row  */

	/*

	   This routine works as follows: 

		allocate temp a vector of size mloc 
		for j=1:n,                                % for all the columns of A
			temp  <- Asave(:,j); 
			normA <- normA + norm(temp,'fro')^2;
			temp  <- temp - Q(:,1:j)*R(1:j,j); % this is a DGEMV operation
			normR <- normR + norm(temp,'fro')^2;
		end
		MPI_Allreduce on normA and normR, the

	   A variant is to perform A-QR with DTRMM but this mean allocating
	   another matrix or size A (eq. Q) and we do not want to do this

	*/
	
	int j;
	double temppp[2];
	double *temp=NULL;

	temp = (double *)malloc(mloc*1*sizeof(double)) ;

	temppp[0] = 0.0e+00;
	temppp[1] = 0.0e+00;
	for ( j = 0; j < n; j++) {
		cblas_dcopy( mloc, &(A[j*mloc]), 1, temp, 1);
		temppp[0] += cblas_ddot(mloc,temp,1,temp,1);
		cblas_dgemv(CblasColMajor, CblasNoTrans, mloc, (j+1), 1,  Q, mloc, &(R[j*n]), 1,
				(-1.0e+00), temp, 1);
		temppp[1] += cblas_ddot(mloc,temp,1,temp,1);
	}
	MPI_Allreduce( MPI_IN_PLACE, &temppp, 2, MPI_DOUBLE, MPI_SUM, mpi_comm);

	free(temp);

	return sqrt(temppp[1]/temppp[0]) ;

}
