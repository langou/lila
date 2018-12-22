#include "project.h"

#define Q(i,j) (Q[ j*lda+i ] )
#define ORTH_TF(i,j) (ORTH_TF[ j*ldatf+i ] )

double check_orthogonality_RFP(int mloc, int n, double *Q, int lda, MPI_Comm mpi_comm){

	double *ORTH_TF = NULL;
	int ldatf;
	int j, my_rank, pool_size;	
	double orthlevel;

	MPI_Comm_size( mpi_comm, &pool_size);
	MPI_Comm_rank( mpi_comm, &my_rank);

	/*
	   We are using RFP format FORM='T', UPLO='L'. 
	   So we store APF = [ U_{22}^T \\ U{11} ; U_{12} ].
	   If n is even, APF is of size (n/2)-by-(n+1).
	   If n is odd,  APF is of size ((n+1)/2)-by-(n).
	   (See picture for more info.)
	*/

	ORTH_TF = (double *)malloc(((n+2)*n)/2*sizeof(double)) ;
	for (j = 0; j < (n*(n+1))/2; j++) ORTH_TF[j]=0xDEADBEEF;

	if (n%2) ldatf = (n+1)/2; else ldatf = n/2;

	if (n%2) { /* this is the case when n is odd */

		/* (1) compute U11 and put in ORTH_TF starting on position ORTH_TF(1,2) */
		cblas_dsyrk( CblasColMajor, CblasUpper, CblasTrans, (n-1)/2, mloc, 1.0e+00, &(Q(0,0)), lda, 0e+00, &(ORTH_TF(0,0)), ldatf); 

		/* (2) compute U11 and put in ORTH_TF starting on position ORTH_TF(1,2) */
		cblas_dgemm( CblasColMajor, CblasTrans, CblasNoTrans,
			(n+1)/2,(n+1)/2,mloc,1e+00,&(Q(0,0)),lda,&(Q(0,(n-1)/2)),lda,0.0e+00,&(ORTH_TF[((n-1)/2)*((n+1)/2)]),ldatf);

		/* (3) compute U22^T and put in ORTH_TF starting on position ORTH_TF(1,1) */
		cblas_dsyrk( CblasColMajor, CblasLower, CblasTrans, (n-1)/2, mloc, 1.0e+00, &(Q(0,(n+1)/2)), lda, 0e+00, &(ORTH_TF(1,0)), ldatf); 


	} else {   /* this is the case when n is even */

		/* (1) compute U11 and put in ORTH_TF starting on position ORTH_TF(1,2) */
		cblas_dsyrk( CblasColMajor, CblasUpper, CblasTrans, n/2, mloc, 1.0e+00, &(Q(0,0)), lda, 0e+00, &(ORTH_TF(0,1)), n/2); 

		/* (2) compute U11 and put in ORTH_TF starting on position ORTH_TF(1,2) */
		cblas_dgemm( CblasColMajor, CblasTrans, CblasNoTrans,
			n/2,n/2,mloc,1e+00,&(Q(0,0)),lda,&(Q(0,n/2)),lda,0.0e+00,&(ORTH_TF[(n/2+1)*(n/2)]),n/2);

		/* (3) compute U22^T and put in ORTH_TF starting on position ORTH_TF(1,1) */
		cblas_dsyrk( CblasColMajor, CblasLower, CblasTrans, n/2, mloc, 1.0e+00, &(Q(0,n/2)), lda, 0e+00, &(ORTH_TF(0,0)), n/2); 

	}

	if (my_rank==0) MPI_Reduce( MPI_IN_PLACE, ORTH_TF, (n*(n+1))/2, MPI_DOUBLE, MPI_SUM, 0, mpi_comm);
	else MPI_Reduce( ORTH_TF, NULL, (n*(n+1))/2, MPI_DOUBLE, MPI_SUM, 0, mpi_comm);

	if (my_rank==0){

		if (n%2) { /* case when n is odd */
			
			for ( j = 0; j < ((n-1)/2); j++) { ORTH_TF[j*((n+1)/2+1)] -= 1.0e+00; ORTH_TF[j*((n+1)/2+1)+1] -= 1.0e+00; }
			ORTH_TF[((n+1)/2)*((n+1)/2)-1] -= 1.0e+00;

		} else {   /* case when n is even */

			ORTH_TF[0] -= 1.0e+00;
			for ( j = 1; j < (n/2); j++) { ORTH_TF[j*(n/2+1)-1] -= 1.0e+00; ORTH_TF[j*(n/2+1)] -= 1.0e+00; }
			ORTH_TF[(n/2-1)*(n/2+1)+n/2] -= 1.0e+00;

		}

	}

	/*
	well there should be soon an LAPACK routine to compute the FROB-norm of an RFP matrix and this might give:

 			orthlevel = dlansf_("F", "T", "L", &n, ORTH_TF, NULL);

	but for the moment we deal with it ourselves
	*/

	if (my_rank==0){
		/* Note that this is quick a dirty: we do not care about possible overflow */

		if (n%2) { /* case when n is odd */
			orthlevel =+ pow(lapack_dlansy( 'F', 'U', (n+1)/2, &(ORTH_TF(0,0)), (n+1)/2, NULL ),2.0e+00);
			orthlevel += pow(lapack_dlansy( 'F', 'L', (n-1)/2, &(ORTH_TF(1,0)), (n+1)/2, NULL ),2.0e+00);
			orthlevel += 2*pow(lapack_dlange( 'F', (n+1)/2, (n-1)/2, &(ORTH_TF[((n+1)/2)*((n+1)/2)]), (n+1)/2, NULL ),2.0e+00);
			orthlevel = sqrt(orthlevel) ;
		}

		else {     /* case when n is even */ 
			orthlevel =+ pow(lapack_dlansy( 'F', 'U', n/2, &(ORTH_TF(0,1)), n/2, NULL ),2.0e+00);
			orthlevel += pow(lapack_dlansy( 'F', 'L', n/2, &(ORTH_TF(0,0)), n/2, NULL ),2.0e+00);
			orthlevel += 2*pow(lapack_dlange( 'F', n/2, n/2, &(ORTH_TF[(n/2+1)*(n/2)]), n/2, NULL ),2.0e+00);
			orthlevel = sqrt(orthlevel) ;
		}
	}


	free(ORTH_TF);
	return orthlevel;

}
