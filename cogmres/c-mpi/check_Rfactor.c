#include "project.h"

int check_Rfactor( double *orthlevel, double *represlevel, int mloc, int n, double *A, int lda, double *R, int ldr, MPI_Comm mpi_comm){

	/*
	  implementation details:
	  =======================
	  this routines allocates: WAY TOO MUCH
	  this routines calls householderA to perform the QR factor of A*R' , could be something else in practice
	*/

	/*
	 * We want to check the backward stability of the R factor, for this we want to solve the orthogonal Procrustes problem
	 *      Find U such that || A - UZ || = min { || A - QR ||, Q has orthonormal columns }
	 *
	 * The way to find U is expensive !!! But well, if you want to check that's the only way I have in mind.
	 *
	 * Basically in Matlab:

	   % input A and R
	   % output Q, the Q factor
	   % (the SVD of A*R' is performed by first a QR facto, then an SVD since A is tall and thin).

	   [m,n]=size(A);
	   Z = A * R';
	   [QZ,RZ] = qr( Z, 0 );
	   [UZ,SZ,VZ] = svd( RZ );
	   Q = QZ(1:m,1:n) * (UZ(1:n,1:n) * VZ(1:n,1:n)');

	   In other words, the problem Q R = A does not determine Q uniquely; and, among all the solutions, we would like a Q that has orthogonal columns

	 */

	double *Z, *UZ, *SZ, *VTZ, *work, *U, *RZ;
	double *Q;
	double ldwork;
	int lwork, info;

	Z = (double *)malloc(mloc*n*sizeof(double)) ;
	RZ  = (double *)malloc(n*n*sizeof(double)) ;
	UZ  = (double *)malloc(n*n*sizeof(double)) ;
	SZ  = (double *)malloc(n*sizeof(double)) ;
	VTZ = (double *)malloc(n*n*sizeof(double)) ;
	U = (double *)malloc(n*n*sizeof(double)) ;
	Q = (double *)malloc(mloc*n*sizeof(double)) ;

	lapack_dlacpy( 'A', mloc, n, A, lda, Z, mloc );
	cblas_dtrmm( CblasColMajor, CblasRight, CblasUpper, CblasTrans, CblasNonUnit,
		mloc, n, 1.0e+00, R, ldr, Z, mloc);

	reducehouseholder_A_v0( mloc, n, Z, mloc, RZ, n, mpi_comm);
	/* I am not sure whether the subdiagonal elements should be zeros or not ... */
	/* So waiting for this to be settled, we set them explicitly to 0.0e+00      */
	{int i,j; for (j=0;j<n;j++) for (i=j+1;i<n;i++) RZ[i+j*n]=0.0e+00; }

	lwork = -1;
	lapack_dgesvd( 'A', 'A', n, n, RZ, ldr, SZ, UZ, n, VTZ, n, &ldwork, lwork, &info );
	lwork = (int) ldwork;
	work = (double *)malloc(lwork*sizeof(double)) ;
	lapack_dgesvd( 'A', 'A', n, n, RZ, ldr, SZ, UZ, n, VTZ, n, work, lwork, &info );
	free(work);

	cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0e+00, UZ, n, VTZ, n, 0.0e+00, U, n);
	cblas_dgemm( CblasColMajor, CblasNoTrans, CblasNoTrans, mloc, n, n, 1.0e+00, Z, mloc, U, n, 0.0e+00, Q, mloc);

	(*orthlevel) = check_orthogonality( mloc, n, Q, mloc, MPI_COMM_WORLD);
	(*represlevel) = check_representativity( mloc, n, A, mloc, Q, mloc, R, n, MPI_COMM_WORLD);

	free(Q);
	free(U);
	free(VTZ);
	free(SZ);
	free(UZ);
	free(RZ);
	free(Z);

	return 0;

}
