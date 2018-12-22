#include "project.h"

int main(int argc, char **argv){
	int i, j, index, m, mloc, my_rank, n, pool_size, seed, verbose;
	double *A=NULL, *Asave=NULL, orthlevel, represlevel, *R=NULL, elapsed_time, orthlevelR, represlevelR;
	char charlilamethod[11];
	int lilamethod;
	int LILA_METHOD_CGS_V0         = 2700;
	int LILA_METHOD_CGS_V1         = 2701;
	int LILA_METHOD_MGS_V0         = 2800;
	int LILA_METHOD_MGS_V1         = 2801;
	int LILA_METHOD_ROWMGS_V0      = 2900;
	int LILA_METHOD_ROWMGS_V1      = 2901;
	int LILA_METHOD_QR2A           = 3000;
	int LILA_METHOD_QR2B           = 3100;
	int LILA_METHOD_QRFA           = 3200;
	int LILA_METHOD_QRFB           = 3300;
	int LILA_METHOD_CHOLQRA_V0     = 3500;
	int LILA_METHOD_CHOLQRA_V1     = 3501;
	int LILA_METHOD_CHOLQRA_V2     = 3502;
	int LILA_METHOD_CHOLQRA_V3     = 3503;
	int LILA_METHOD_CHOLQRB_V0     = 3600;
	int LILA_METHOD_CHOLQRB_V1     = 3601;
	int LILA_METHOD_CHOLQRB_V2     = 3602;
	int LILA_METHOD_CHOLQRB_V3     = 3603;
	int LILA_METHOD_RHHA_V0        = 3700;
	int LILA_METHOD_RHHA_V1        = 3701;
	int LILA_METHOD_RHHA_V2        = 3702;
	int LILA_METHOD_RHHB_V0        = 3800;
	int LILA_METHOD_RHHB_V1        = 3801;
	int LILA_METHOD_RHHB_V2        = 3802;
	int LILA_METHOD_RHHB_V3        = 3803;
	int (*lila_method_orthogonalization)( int, int, double *, int , double *, int , MPI_Comm) = NULL;

        MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &pool_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	
	m = 100; n = 10; verbose = 1; lilamethod = LILA_METHOD_CHOLQRB_V0;

	for( i = 1; i < argc; i++ ) {
		if( strcmp( argv[i], "-m" ) == 0 ) { m  = atoi(argv[i+1]); i++; }
		if( strcmp( argv[i], "-n" ) == 0 ) { n  = atoi(argv[i+1]); i++; }
		if( strcmp( argv[i], "-verbose" ) == 0 ) { verbose  = atoi(argv[i+1]); i++; }
		if( strcmp( argv[i], "-cgs_v0" ) == 0 ) lilamethod = LILA_METHOD_CGS_V0;
		if( strcmp( argv[i], "-cgs_v1" ) == 0 ) lilamethod = LILA_METHOD_CGS_V1;
		if( strcmp( argv[i], "-mgs_v0" ) == 0 ) lilamethod = LILA_METHOD_MGS_V0;
		if( strcmp( argv[i], "-mgs_v1" ) == 0 ) lilamethod = LILA_METHOD_MGS_V1;
		if( strcmp( argv[i], "-rowmgs_v0" ) == 0 ) lilamethod = LILA_METHOD_ROWMGS_V0;
		if( strcmp( argv[i], "-rowmgs_v1" ) == 0 ) lilamethod = LILA_METHOD_ROWMGS_V1;
		if( strcmp( argv[i], "-qr2A" ) == 0 ) lilamethod = LILA_METHOD_QR2A;
		if( strcmp( argv[i], "-qr2B" ) == 0 ) lilamethod = LILA_METHOD_QR2B;
		if( strcmp( argv[i], "-qrfA" ) == 0 ) lilamethod = LILA_METHOD_QRFA;
		if( strcmp( argv[i], "-qrfB" ) == 0 ) lilamethod = LILA_METHOD_QRFB;
		if( strcmp( argv[i], "-cholqrA_v0" ) == 0 ) lilamethod = LILA_METHOD_CHOLQRA_V0;
		if( strcmp( argv[i], "-cholqrA_v1" ) == 0 ) lilamethod = LILA_METHOD_CHOLQRA_V1;
		if( strcmp( argv[i], "-cholqrA_v2" ) == 0 ) lilamethod = LILA_METHOD_CHOLQRA_V2;
		if( strcmp( argv[i], "-cholqrA_v3" ) == 0 ) lilamethod = LILA_METHOD_CHOLQRA_V3;
		if( strcmp( argv[i], "-cholqrB_v0" ) == 0 ) lilamethod = LILA_METHOD_CHOLQRB_V0;
		if( strcmp( argv[i], "-cholqrB_v1" ) == 0 ) lilamethod = LILA_METHOD_CHOLQRB_V1;
		if( strcmp( argv[i], "-cholqrB_v2" ) == 0 ) lilamethod = LILA_METHOD_CHOLQRB_V2;
		if( strcmp( argv[i], "-cholqrB_v3" ) == 0 ) lilamethod = LILA_METHOD_CHOLQRB_V3;
		if( strcmp( argv[i], "-rhhA_v0" ) == 0 ) lilamethod = LILA_METHOD_RHHA_V0;
		if( strcmp( argv[i], "-rhhA_v1" ) == 0 ) lilamethod = LILA_METHOD_RHHA_V1;
		if( strcmp( argv[i], "-rhhA_v2" ) == 0 ) lilamethod = LILA_METHOD_RHHA_V2;
		if( strcmp( argv[i], "-rhhB_v0" ) == 0 ) lilamethod = LILA_METHOD_RHHB_V0;
		if( strcmp( argv[i], "-rhhB_v1" ) == 0 ) lilamethod = LILA_METHOD_RHHB_V1;
		if( strcmp( argv[i], "-rhhB_v2" ) == 0 ) lilamethod = LILA_METHOD_RHHB_V2;
		if( strcmp( argv[i], "-rhhB_v3" ) == 0 ) lilamethod = LILA_METHOD_RHHB_V3;
	}

	if (lilamethod == LILA_METHOD_CGS_V0)         strcpy( charlilamethod, "CGS_V0     ");
	if (lilamethod == LILA_METHOD_CGS_V1)         strcpy( charlilamethod, "CGS_V1     ");
	if (lilamethod == LILA_METHOD_MGS_V0)         strcpy( charlilamethod, "MGS_V0     ");
	if (lilamethod == LILA_METHOD_MGS_V1)         strcpy( charlilamethod, "MGS_V1     ");
	if (lilamethod == LILA_METHOD_ROWMGS_V0)      strcpy( charlilamethod, "ROWMGS_V0  ");
	if (lilamethod == LILA_METHOD_ROWMGS_V1)      strcpy( charlilamethod, "ROWMGS_V1  ");
	if (lilamethod == LILA_METHOD_QR2A)           strcpy( charlilamethod, "QR2A       ");
	if (lilamethod == LILA_METHOD_QR2B)           strcpy( charlilamethod, "QR2B       ");
	if (lilamethod == LILA_METHOD_QRFA)           strcpy( charlilamethod, "QRFA       ");
	if (lilamethod == LILA_METHOD_QRFB)           strcpy( charlilamethod, "QRFB       ");
	if (lilamethod == LILA_METHOD_CHOLQRA_V0)     strcpy( charlilamethod, "CHOLQRA_V0 ");
	if (lilamethod == LILA_METHOD_CHOLQRA_V1)     strcpy( charlilamethod, "CHOLQRA_V1 ");
	if (lilamethod == LILA_METHOD_CHOLQRA_V2)     strcpy( charlilamethod, "CHOLQRA_V2 ");
	if (lilamethod == LILA_METHOD_CHOLQRA_V3)     strcpy( charlilamethod, "CHOLQRA_V3 ");
	if (lilamethod == LILA_METHOD_CHOLQRB_V0)     strcpy( charlilamethod, "CHOLQRB_V0 ");
	if (lilamethod == LILA_METHOD_CHOLQRB_V1)     strcpy( charlilamethod, "CHOLQRB_V1 ");
	if (lilamethod == LILA_METHOD_CHOLQRB_V2)     strcpy( charlilamethod, "CHOLQRB_V2 ");
	if (lilamethod == LILA_METHOD_CHOLQRB_V3)     strcpy( charlilamethod, "CHOLQRB_V3 ");
	if (lilamethod == LILA_METHOD_RHHA_V0)        strcpy( charlilamethod, "RHHA_V0    ");
	if (lilamethod == LILA_METHOD_RHHA_V1)        strcpy( charlilamethod, "RHHA_V1    ");
	if (lilamethod == LILA_METHOD_RHHA_V2)        strcpy( charlilamethod, "RHHA_V2    ");
	if (lilamethod == LILA_METHOD_RHHB_V0)        strcpy( charlilamethod, "RHHB_V0    ");
	if (lilamethod == LILA_METHOD_RHHB_V1)        strcpy( charlilamethod, "RHHB_V1    ");
	if (lilamethod == LILA_METHOD_RHHB_V2)        strcpy( charlilamethod, "RHHB_V2    ");
	if (lilamethod == LILA_METHOD_RHHB_V3)        strcpy( charlilamethod, "RHHB_V3    ");

	if (lilamethod == LILA_METHOD_CGS_V0)         lila_method_orthogonalization = cgs_v0;
	if (lilamethod == LILA_METHOD_CGS_V1)         lila_method_orthogonalization = cgs_v1;
	if (lilamethod == LILA_METHOD_MGS_V0)         lila_method_orthogonalization = mgs_v0;
	if (lilamethod == LILA_METHOD_MGS_V1)         lila_method_orthogonalization = mgs_v1;
	if (lilamethod == LILA_METHOD_ROWMGS_V0)      lila_method_orthogonalization = rowmgs_v0;
	if (lilamethod == LILA_METHOD_ROWMGS_V1)      lila_method_orthogonalization = rowmgs_v1;
    //	if (lilamethod == LILA_METHOD_QR2A)           lila_method_orthogonalization = scalapackqr2_A;
    //	if (lilamethod == LILA_METHOD_QR2B)           lila_method_orthogonalization = scalapackqr2_B;
    //	if (lilamethod == LILA_METHOD_QRFA)           lila_method_orthogonalization = scalapackqrf_A;
    //	if (lilamethod == LILA_METHOD_QRFB)           lila_method_orthogonalization = scalapackqr2_B;
	if (lilamethod == LILA_METHOD_CHOLQRA_V0)     lila_method_orthogonalization = choleskyqr_A_v0;
	if (lilamethod == LILA_METHOD_CHOLQRA_V1)     lila_method_orthogonalization = choleskyqr_A_v1;
	if (lilamethod == LILA_METHOD_CHOLQRA_V2)     lila_method_orthogonalization = choleskyqr_A_v2;
	if (lilamethod == LILA_METHOD_CHOLQRA_V3)     lila_method_orthogonalization = choleskyqr_A_v3;
	if (lilamethod == LILA_METHOD_CHOLQRB_V0)     lila_method_orthogonalization = choleskyqr_B_v0;
	if (lilamethod == LILA_METHOD_CHOLQRB_V1)     lila_method_orthogonalization = choleskyqr_B_v1;
	if (lilamethod == LILA_METHOD_CHOLQRB_V2)     lila_method_orthogonalization = choleskyqr_B_v2;
	if (lilamethod == LILA_METHOD_CHOLQRB_V3)     lila_method_orthogonalization = choleskyqr_B_v3;
	if (lilamethod == LILA_METHOD_RHHA_V0)        lila_method_orthogonalization = reducehouseholder_A_v0;
	if (lilamethod == LILA_METHOD_RHHA_V1)        lila_method_orthogonalization = reducehouseholder_A_v1;
	if (lilamethod == LILA_METHOD_RHHA_V2)        lila_method_orthogonalization = reducehouseholder_A_v2;
	if (lilamethod == LILA_METHOD_RHHB_V0)        lila_method_orthogonalization = reducehouseholder_B_v0;
	if (lilamethod == LILA_METHOD_RHHB_V1)        lila_method_orthogonalization = reducehouseholder_B_v1;
	if (lilamethod == LILA_METHOD_RHHB_V2)        lila_method_orthogonalization = reducehouseholder_B_v2;
	if (lilamethod == LILA_METHOD_RHHB_V3)        lila_method_orthogonalization = reducehouseholder_B_v3;

	seed = my_rank*m*n; srand(seed);

	mloc = m / pool_size; if ( my_rank < m - pool_size*mloc ) mloc++;

	A = (double *)malloc(mloc*n*sizeof(double)) ;
	Asave = (double *)malloc(mloc*n*sizeof(double)) ;
	R = (double *)malloc(n*n*sizeof(double)) ;

	index = 0;
	for (i = 0; i < mloc; i++) {
		for (j = 0; j < n; j++) {
			A[index] = ((double) rand()) / ((double) RAND_MAX) - 0.5 ;
			Asave[index] = A[index] ;
			index++;	
		}
	}

	/* Initialization of R, not mandatorry a priori, so we in purpose put weird values */
	for (i = 0; i < n*n; i++)  R[i] = 0xDEADBEEF;

	/*
	{
		FILE *file=NULL;
		int myrank;
		int ir, jr;

		MPI_Comm_rank( MPI_COMM_WORLD, &myrank); 
		if ( myrank == 0) { 
			file = fopen ("MATLAB_A1.dat","w");
			for (ir=0;ir<mloc;ir++) {
				for (jr=0;jr<n;jr++){
					fprintf( file, " %+16.14e", A[jr*mloc + ir]);
				}
				fprintf( file, "\n");
			}
			fclose( file );
		}
		if ( myrank == 1) { 
			file = fopen ("MATLAB_A2.dat","w");
			for (ir=0;ir<mloc;ir++) {
				for (jr=0;jr<n;jr++){
					fprintf( file, " %+16.14e", A[jr*mloc + ir]);
				}
				fprintf( file, "\n");
			}
			fclose( file );
		}
		if ( myrank == 2) { 
			file = fopen ("MATLAB_A3.dat","w");
			for (ir=0;ir<mloc;ir++) {
				for (jr=0;jr<n;jr++){
					fprintf( file, " %+16.14e", A[jr*mloc + ir]);
				}
				fprintf( file, "\n");
			}
			fclose( file );
		}
		if ( myrank == 3) { 
			file = fopen ("MATLAB_A4.dat","w");
			for (ir=0;ir<mloc;ir++) {
				for (jr=0;jr<n;jr++){
					fprintf( file, " %+16.14e", A[jr*mloc + ir]);
				}
				fprintf( file, "\n");
			}
			fclose( file );
		}
	}
	*/

	/* Call the orthogonalization scheme */
	MPI_Barrier(MPI_COMM_WORLD);
	elapsed_time = -MPI_Wtime();
	(*lila_method_orthogonalization)( mloc, n, A, mloc, R, n, MPI_COMM_WORLD);
	elapsed_time += MPI_Wtime();

	/* Checks and Prints */

	if ( (lilamethod == LILA_METHOD_QR2B)
	  || (lilamethod == LILA_METHOD_QRFB)
	  || (lilamethod == LILA_METHOD_CHOLQRB_V0 )
	  || (lilamethod == LILA_METHOD_CHOLQRB_V1 )
	  || (lilamethod == LILA_METHOD_CHOLQRB_V2 )
	  || (lilamethod == LILA_METHOD_CHOLQRB_V3 )
	  || (lilamethod == LILA_METHOD_RHHB_V0 ) 
	  || (lilamethod == LILA_METHOD_RHHB_V1 ) 
	  || (lilamethod == LILA_METHOD_RHHB_V2 )
	  || (lilamethod == LILA_METHOD_RHHB_V3 ) ) {

		{
			int lda=mloc,ldr=n;
			for (j = 0; j < n; j++) for (i = 0; i < mloc; i++) A[lda*j+i]=Asave[lda*j+i];
			cblas_dtrsm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, mloc, n, 1.0e+00, R, ldr, A, lda );
			orthlevel = check_orthogonality( mloc, n, A, lda, MPI_COMM_WORLD);
			represlevel = check_representativity( mloc, n, Asave, mloc, A, lda, R, ldr, MPI_COMM_WORLD);
		}
		{
			int lda=mloc,ldr=n;
			check_Rfactor( &orthlevelR, &represlevelR, mloc, n, Asave, lda, R, ldr, MPI_COMM_WORLD);
		}
		if (my_rank==0) {
			if (verbose==0) { }
			else if (verbose==2){
				printf("\n%s :: #proc = %d, m = %d, n = %d\n",charlilamethod, pool_size, m, n);
				printf("%s :: elapsed time          =  %10.3f sec\n",charlilamethod,elapsed_time);
				printf("%s :: Flops rate (total)    =  %10.3f MFlops/sec\n",charlilamethod,
					2.00e+00*((double) m)*((double) n)*((double) n)*1.00e-6/elapsed_time);
				printf("%s :: Flops rate / proc     =  %10.3f MFlops/sec/proc\n",charlilamethod,
					2.00e+00*((double) m)*((double) n)*((double) n)*1.00e-6/elapsed_time/pool_size);
				printf("%s :: orthlevel             =  %e\n",charlilamethod,orthlevel);
				printf("%s :: represlevel           =  %e\n",charlilamethod,represlevel);
				printf("%s :: orthlevelR            =  %e\n",charlilamethod,orthlevelR);
				printf("%s :: represlevelR          =  %e\n",charlilamethod,represlevelR);
				printf("\n");
			}
			else {
				printf("\t%4d\t%4d\t%12d\t%6d\t%10.3f\t%10.3f\t%10.3f\t%6.0e\t%6.0e\t %6.0e\t %6.0e\n",
					lilamethod,
					pool_size, m, n,
					elapsed_time,
					2.00e+00*((double) m)*((double) n)*((double) n)*1.00e-6/elapsed_time,
					2.00e+00*((double) m)*((double) n)*((double) n)*1.00e-6/elapsed_time/pool_size,
					orthlevel, represlevel,
					orthlevelR, represlevelR );
			}
		}
	}
	else {
		orthlevel = check_orthogonality( mloc, n, A, mloc, MPI_COMM_WORLD);
		represlevel = check_representativity( mloc, n, Asave, mloc, A, mloc, R, n, MPI_COMM_WORLD);
		if (my_rank==0) {
			if (verbose==0) { }
			else if (verbose==2){
				printf("\n%s :: #proc = %d, m = %d, n = %d\n",charlilamethod, pool_size, m, n);
				printf("%s :: elapsed time          =  %10.3f sec\n",charlilamethod,elapsed_time);
				printf("%s :: Flops rate (total)    =  %10.3f MFlops/sec\n",charlilamethod,
					2.00e+00*((double) m)*((double) n)*((double) n)*1.00e-6/elapsed_time);
				printf("%s :: Flops rate / proc     =  %10.3f MFlops/sec/proc\n",charlilamethod,
					2.00e+00*((double) m)*((double) n)*((double) n)*1.00e-6/elapsed_time/pool_size);
				printf("%s :: orthlevel             =  %e\n",charlilamethod,orthlevel);
				printf("%s :: represlevel           =  %e\n",charlilamethod,represlevel);
				printf("\n");
			}
			else {
				printf("\t%4d\t%4d\t%12d\t%6d\t%10.3f\t%10.3f\t%10.3f\t%6.0e\t%6.0e\t NaN\t NaN\n",
					lilamethod,
					pool_size, m, n,
					elapsed_time,
					2.00e+00*((double) m)*((double) n)*((double) n)*1.00e-6/elapsed_time,
					2.00e+00*((double) m)*((double) n)*((double) n)*1.00e-6/elapsed_time/pool_size,
					orthlevel, represlevel);
			}
		}
	}

	free(R);
	free(Asave);
	free(A);

        MPI_Finalize();
	exit(0);

}
