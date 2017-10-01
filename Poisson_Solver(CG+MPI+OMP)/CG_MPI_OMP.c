#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "mpi.h"
#include "omp.h"

#define PI 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651328
#define StopEps 1e-6							// Stop Condition.
#define MaxIter 10000							// Maximum #Iteration. 
#define DefineNumOfGrid 256						// To use #Grid.
#define DebugFlag 1								// 0: off, 1: on.

double *AnalyticalSol_for_PurePoisson(int NumOfGrid);
double *ForceTerm_for_PurePoisson(int NumOfGrid);
double *CG_for_PurePoisson(int argc, char * argv[], double *x0, double *Fvalue, int NumOfGrid);

////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////     Function	(Run)	////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
double *AnalyticalSol_for_PurePoisson(int NumOfGrid) {
	int i, j;											// Iteration Index.
	double dxdygrid = (NumOfGrid - 1);					// NumOfGrid for loop.
	double *ASol;										// Storage of Analytical Solution.

														// Initialize array using malloc fucntion
	ASol = (double*)malloc(sizeof(double)*NumOfGrid*NumOfGrid);

	// If you want to input a Analytical Soluton, you modify the equation below :
	for (i = 0; i < NumOfGrid; i++)
		for (j = 0; j < NumOfGrid; j++)
			ASol[j + i*NumOfGrid] = -0.5*(sin(PI*i / dxdygrid)*cos(PI*j / dxdygrid) / (PI*PI));

	return ASol;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
double *ForceTerm_for_PurePoisson(int NumOfGrid) {
	int i, j;											// Iteration Index.
	double dxdygrid = (NumOfGrid - 1);					// NumOfGrid for loop.
	double *Fvalue;									// Storage of Analytical Solution.

													// Initialize array using malloc fucntion
	Fvalue = (double*)malloc(sizeof(double)*NumOfGrid*NumOfGrid);

	// If you want to input a Force Term, you modify the equation below :
	for (i = 0; i < NumOfGrid; i++)
		for (j = 0; j < NumOfGrid; j++)
			Fvalue[j + i*NumOfGrid] = (sin(PI*i / dxdygrid)*cos(PI*j / dxdygrid)) / ((dxdygrid + 1)*(dxdygrid + 1));

	return Fvalue;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
double *CG_for_PurePoisson(int argc, char * argv[], double *x0, double *Fvalue, int NumOfGrid) {

	MPI_Init(&argc, &argv);

	int nrank;
	int nproc;

	MPI_Comm_rank(MPI_COMM_WORLD, &nrank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Barrier(MPI_COMM_WORLD);

	if (((double)sqrt(nproc) - (int)sqrt(nproc)) != 0) {
		printf("Error - Do not decompose Matrix to squres. \n");
		exit(-1);
	}

	int GridMPI, root = 0;
	int i, j, k;													 // Iteration Index.
	int ITER = 0, loop, MPI_ITER;									 // #Iteration.
	double ClockTime;												 // Consumed Time.
	double dxdygrid = (NumOfGrid - 1);								 // NumOfGrid for loop.
	double N2Err = 0;												 // (Sub) Stop Condition.
	double SumOfN2Err = 1;											 // Initial Stop Condition for (While).
	double *x, *xMPI;												 // Storage of Updated Matrix 
	double *new_x, *new_xMPI;										 // Storage of Updated Matrix 
	double *r0, *r0MPI, *rM0MPI;									 // Storage of OLD Residual Matrix.
	double *p0, *p0MPI, *pM0MPI;									 // Storage of OLD Direction Matrix. 
	double *FvalueVector, *FvalueVectorMPI;							 // Force Term Vector.
	double *PoissonMultipleVector, *PoissonMultipleVectorMPI;		 // Storage of Updated Matrix.
	double Alpha, Beta;												 // CG Parameter.
	double UpperValue = 0, LowerValue = 0;							 // Constant Storage for CG Parameter.
	double UpperValueLocal, LowerValueLocal, SumOfN2ErrLocal, t0, t1;// MPI Reduce.
	double *nnz;													 // CRS Variable.
	int *row_ptr, *col_ind;											 // CRS Variable.
	int CRS_iter = 3, tid;
	int rowptriter;
	
	x = (double*)malloc(sizeof(double)*NumOfGrid*NumOfGrid);
	new_x = (double*)malloc(sizeof(double)*NumOfGrid*NumOfGrid);
	MPI_Request request, request2;
	MPI_Status status;
	// CRS array
	nnz = (double*)malloc(sizeof(double)*NumOfGrid*NumOfGrid * 5);
	row_ptr = (int*)malloc(sizeof(int)*NumOfGrid*NumOfGrid);
	col_ind = (int*)malloc(sizeof(int)*NumOfGrid*NumOfGrid * 5);
	
	r0 = (double*)malloc(sizeof(double)*NumOfGrid*NumOfGrid);
	p0 = (double*)malloc(sizeof(double)*NumOfGrid*NumOfGrid);
		
	
	nnz[0] = -4.0;
	nnz[1] = 2.0;
	nnz[2] = 1.0;
	row_ptr[0] = 0;
	col_ind[0] = 0;
	col_ind[1] = 1;
	col_ind[2] = NumOfGrid;

	for (i = 2; i <= NumOfGrid*NumOfGrid - 1; i++) {
		row_ptr[i - 1] = CRS_iter;
		for (j = 1; j <= (NumOfGrid*NumOfGrid); j++) {
			if (i == j) {
				nnz[CRS_iter] = -4;
				col_ind[CRS_iter] = j - 1;
				CRS_iter += 1;
			}
			else if (((i + 1) == j) && ((i % (NumOfGrid)) != 0)) {
				if (((i - 1) % (NumOfGrid)) == 0)
					nnz[CRS_iter] = 2;
				else
					nnz[CRS_iter] = 1;
				col_ind[CRS_iter] = j - 1;
				CRS_iter += 1;
			}
			else if ((i == (j + 1)) && ((j % (NumOfGrid)) != 0)) {
				if ((i % (NumOfGrid)) == 0)
					nnz[CRS_iter] = 2;
				else
					nnz[CRS_iter] = 1;
				col_ind[CRS_iter] = j - 1;
				CRS_iter += 1;
			}
			else if ((i + NumOfGrid) == j) {
				nnz[CRS_iter] = 1;
				col_ind[CRS_iter] = j - 1;
				CRS_iter += 1;
			}
			else if (i == (j + NumOfGrid)) {
				nnz[CRS_iter] = 1;
				col_ind[CRS_iter] = j - 1;
				CRS_iter += 1;
			}
		}
	}
	row_ptr[NumOfGrid*NumOfGrid - 1] = CRS_iter;
	row_ptr[NumOfGrid*NumOfGrid] = CRS_iter + 3;
	nnz[CRS_iter] = 1;
	nnz[CRS_iter + 1] = 2;
	nnz[CRS_iter + 2] = -4;
	col_ind[CRS_iter] = NumOfGrid*NumOfGrid - NumOfGrid - 1;
	col_ind[CRS_iter + 1] = NumOfGrid*NumOfGrid - 2;
	col_ind[CRS_iter + 2] = NumOfGrid*NumOfGrid - 1;



	int CRS_Size = CRS_iter + 3, row_num = 0;
	int MPICRS_ColStart, MPICRS_ColEnd, MPICRS_RowStart, MPICRS_RowEnd, *PARTrow_ptr, *PARTcol_ind, *row_ind, *PARTrow_ind;
	double *PARTnnz;

	GridMPI = NumOfGrid*NumOfGrid / sqrt(nproc);

	PARTnnz = (double*)malloc(sizeof(double) * GridMPI * 5);
	PARTrow_ptr = (int*)malloc(sizeof(int) * GridMPI);
	PARTcol_ind = (int*)malloc(sizeof(int) * GridMPI * 5);
	PARTrow_ind = (int*)malloc(sizeof(int) * GridMPI * 5);
	row_ind = (int*)malloc(sizeof(int) * NumOfGrid*NumOfGrid * 5);

	for (i = 0; i < (NumOfGrid*NumOfGrid); i++) {
		for (j = 0; j < (row_ptr[i + 1] - row_ptr[i]); j++) {
			row_ind[ITER] = row_num;
			ITER = ITER + 1;
		}
		row_num = row_num + 1;
	}


	if (nrank < sqrt(nproc)) {
		MPICRS_ColStart = GridMPI * nrank;
		MPICRS_ColEnd = (GridMPI)* nrank + GridMPI;
		MPICRS_RowStart = 0;
		MPICRS_RowEnd = GridMPI;
		int Destination = nrank + sqrt(nproc);
		MPI_Isend(&MPICRS_ColStart, 1, MPI_INT, Destination, 101, MPI_COMM_WORLD, &request);
		MPI_Isend(&MPICRS_ColEnd, 1, MPI_INT, Destination, 102, MPI_COMM_WORLD, &request2);
		
	}
	else if (nrank >= sqrt(nproc)) {
		int RecvSource = nrank - sqrt(nproc);
		int RowPlus = 1;
		MPI_Irecv(&MPICRS_ColStart, 1, MPI_INT, RecvSource, 101, MPI_COMM_WORLD, &request);
		MPI_Irecv(&MPICRS_ColEnd, 1, MPI_INT, RecvSource, 102, MPI_COMM_WORLD, &request2);
		MPICRS_RowStart = GridMPI * (int)sqrt(nrank);
		MPICRS_RowEnd = GridMPI * ((int)sqrt(nrank) + 1);
		MPI_Wait(&request, &status);
		MPI_Wait(&request2, &status);
	}

	CRS_iter = 0;
	for (i = 0; i < CRS_Size; i++)
		if ((col_ind[i] >= MPICRS_ColStart) && (col_ind[i] < MPICRS_ColEnd) &&
			(row_ind[i] >= MPICRS_RowStart) && (row_ind[i] < MPICRS_RowEnd)) {
			PARTcol_ind[CRS_iter] = col_ind[i] - MPICRS_ColStart;
			PARTrow_ind[CRS_iter] = row_ind[i] - MPICRS_RowStart;
			PARTnnz[CRS_iter] = nnz[i];
			CRS_iter = CRS_iter + 1;
		}

	// Initialize array using malloc fucntion (MPI)
	new_xMPI = (double*)malloc(sizeof(double)*GridMPI);
	xMPI = (double*)malloc(sizeof(double)*GridMPI);
	r0MPI = (double*)malloc(sizeof(double)*GridMPI);
	rM0MPI = (double*)malloc(sizeof(double)*GridMPI);
	p0MPI = (double*)malloc(sizeof(double)*GridMPI);
	pM0MPI = (double*)malloc(sizeof(double)*GridMPI);
	PoissonMultipleVectorMPI = (double*)malloc(sizeof(double)*GridMPI);


	for (i = 0; i < GridMPI; i++) {
		xMPI[i] = x0[i + MPICRS_RowStart];
		r0MPI[i] = Fvalue[i + MPICRS_ColStart];
		rM0MPI[i] = Fvalue[i + MPICRS_RowStart];
		p0MPI[i] = r0MPI[i];
		pM0MPI[i] = rM0MPI[i];
	}
	
	
	int RowINDEX = 0, tempValue, t, RowRank = 0, ColRank, RowStep = 0;
	double RightValue = 0, LeftValue = 0, RowSum = 0, iSend = 0, iRecv = 0;
	
	MPI_Barrier(MPI_COMM_WORLD);
	t0 = MPI_Wtime();
	for (loop = 0; loop < MaxIter; loop++) {
		#pragma omp parallel for private(i) num_threads(4)
		for (i = 0; i < NumOfGrid * NumOfGrid; i++)
			new_x[i] = x[i];
		
		UpperValueLocal = 0.0;
		#pragma omp parallel for private(i) reduction(+:UpperValueLocal) num_threads(4)
		for (i = 0; i < GridMPI; i++)
			UpperValueLocal += (pM0MPI[i] * rM0MPI[i]);
		
		UpperValue = 0.0;
		if (nrank == 3)
			MPI_Isend(&UpperValueLocal, 1, MPI_DOUBLE, 0, 103, MPI_COMM_WORLD, &request);
		if (nrank == 0) {
			UpperValue = UpperValueLocal;
			MPI_Irecv(&UpperValueLocal, 1, MPI_DOUBLE, 3, 103, MPI_COMM_WORLD, &request);
			//MPI_Wait(&request, &status);
			UpperValue += UpperValueLocal;}
		MPI_Bcast(&UpperValue, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		#pragma omp parallel for private(i) num_threads(4)
		for (i = 0; i < GridMPI; i++)
			PoissonMultipleVectorMPI[i] = 0.0;

		ITER = 0;
		for (RowRank = 0; RowRank < GridMPI; RowRank++)
		{
			LeftValue = 0.0;
			#pragma omp parallel for private(i) reduction(+:LeftValue) num_threads(4)
			for (i = 0; i < CRS_iter; i++)
				if (PARTrow_ind[i] == RowRank)
					LeftValue += PARTnnz[i] * p0MPI[PARTcol_ind[i]];
					
			for (ColRank = 0; ColRank < sqrt(nproc) - 1; ColRank++) {
				for (RowStep = 0; RowStep <= sqrt(nproc); RowStep += sqrt(nproc)) {
					if (ColRank + RowStep == nrank)
						MPI_Isend(&LeftValue, 1, MPI_DOUBLE, ColRank + RowStep + 1, 104, MPI_COMM_WORLD, &request);
					else if ((ColRank + RowStep + 1) == nrank) {
						RightValue = LeftValue;
						MPI_Irecv(&LeftValue, 1, MPI_DOUBLE, ColRank + RowStep, 104, MPI_COMM_WORLD, &request);
						MPI_Wait(&request, &status);
						LeftValue = RightValue + LeftValue;
					}
				}
			}
			

			for (ColRank = sqrt(nproc) - 1; ColRank > 0; ColRank--)
				for (RowStep = 0; RowStep <= sqrt(nproc); RowStep += sqrt(nproc)) {
					if (ColRank + RowStep == nrank)
						MPI_Isend(&LeftValue, 1, MPI_DOUBLE, ColRank + RowStep - 1, 105, MPI_COMM_WORLD, &request);
					else if ((ColRank + RowStep - 1) == nrank){
						MPI_Irecv(&LeftValue, 1, MPI_DOUBLE, ColRank + RowStep, 105, MPI_COMM_WORLD, &request);
						MPI_Wait(&request, &status);
					}
				}


			PoissonMultipleVectorMPI[ITER] = LeftValue;
			ITER += 1;
		}

	
		LowerValueLocal = 0.0;
		#pragma omp parallel for private(i) reduction(+:LowerValueLocal) num_threads(4)
		for (i = 0; i < GridMPI; i++)
			LowerValueLocal += pM0MPI[i] * PoissonMultipleVectorMPI[i];
		
		if (nrank == 3)
			MPI_Isend(&LowerValueLocal, 1, MPI_DOUBLE, 0, 106, MPI_COMM_WORLD, &request);


		LowerValue = 0.0;
		if (nrank == 0) {
			LowerValue = LowerValueLocal;
			MPI_Irecv(&LowerValueLocal, 1, MPI_DOUBLE, 3, 106, MPI_COMM_WORLD, &request);
			//MPI_Wait(&request, &status);
			LowerValue += LowerValueLocal;}

		MPI_Bcast(&LowerValue, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		
		Alpha = UpperValue / LowerValue;
		

		#pragma omp parallel for private(i) num_threads(4)
		for (i = 0; i < GridMPI; i++) {
			xMPI[i] = xMPI[i] + Alpha * pM0MPI[i];
			rM0MPI[i] = rM0MPI[i] - Alpha * PoissonMultipleVectorMPI[i];}

		UpperValueLocal = 0.0;
		#pragma omp parallel for private(i) reduction(+:UpperValueLocal) num_threads(4)
		for (i = 0; i < GridMPI; i++)
			UpperValueLocal += (rM0MPI[i] * PoissonMultipleVectorMPI[i]);
		
		if (nrank == 3)
			MPI_Isend(&UpperValueLocal, 1, MPI_DOUBLE, 0, 107, MPI_COMM_WORLD, &request);
		
		UpperValue = 0.0;
		if (nrank == 0) {
			UpperValue = UpperValueLocal;
			MPI_Irecv(&UpperValueLocal, 1, MPI_DOUBLE, 3, 107, MPI_COMM_WORLD, &request);
			MPI_Wait(&request, &status);
			UpperValue += UpperValueLocal;}

		MPI_Bcast(&UpperValue, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);


		Beta = -UpperValue / LowerValue;

		#pragma omp parallel for private(i) num_threads(4)
		for (i = 0; i < GridMPI; i++)
			pM0MPI[i] = rM0MPI[i] - Beta * pM0MPI[i];

		if (nrank == 0) {
			#pragma omp parallel for private(i) num_threads(4)
			for (i = 0; i < (NumOfGrid*NumOfGrid); i++)
				x[i] = 0;
			#pragma omp parallel for private(i) num_threads(4)
			for (i = 0; i < GridMPI; i++) 
				x[i] = xMPI[i];}

		if (nrank == 3)
			MPI_Isend(&xMPI[0], GridMPI, MPI_DOUBLE, 0, 108, MPI_COMM_WORLD, &request);
		else if (nrank == 0) {
			MPI_Irecv(&xMPI[0], GridMPI, MPI_DOUBLE, 3, 108, MPI_COMM_WORLD, &request);
			MPI_Wait(&request, &status);
			#pragma omp parallel for private(i) num_threads(4)
			for (i = 0; i < GridMPI; i++) 
				x[i + GridMPI] = xMPI[i];}

		MPI_Bcast(x, (NumOfGrid*NumOfGrid), MPI_DOUBLE, 0, MPI_COMM_WORLD);

		#pragma omp parallel for private(i) num_threads(4)
		for (i = 0; i < GridMPI; i++)
			xMPI[i] = x[i + MPICRS_RowStart];

		SumOfN2Err = 0;
		#pragma omp parallel for private(i) reduction(+:SumOfN2Err) num_threads(4)
		for (i = 0; i < (NumOfGrid*NumOfGrid); i++) {
			N2Err = (new_x[i] - x[i])*(new_x[i] - x[i]);
			SumOfN2Err += N2Err;}

		if (nrank == 0) {
			for (i = 0; i < (NumOfGrid*NumOfGrid); i++) {
				r0[i] = 0;
				p0[i] = 0;}
			for (i = 0; i < GridMPI; i++) {
				r0[i] = rM0MPI[i];
				rM0MPI[i] = 0;
				p0[i] = pM0MPI[i];
				pM0MPI[i] = 0;}	}
		
		if (nrank == 3) {
			MPI_Isend(&rM0MPI[0], GridMPI, MPI_DOUBLE, 0, 507, MPI_COMM_WORLD, &request);
			MPI_Isend(&pM0MPI[0], GridMPI, MPI_DOUBLE, 0, 508, MPI_COMM_WORLD, &request2);
		}
		else if (nrank == 0) {
			MPI_Irecv(&rM0MPI[0], GridMPI, MPI_DOUBLE, 3, 507, MPI_COMM_WORLD, &request);
			MPI_Irecv(&pM0MPI[0], GridMPI, MPI_DOUBLE, 3, 508, MPI_COMM_WORLD, &request2);
			MPI_Wait(&request, &status);
			MPI_Wait(&request2, &status);
			#pragma omp parallel for private(i) num_threads(4)
			for (i = 0; i < GridMPI; i++){
				r0[i + GridMPI] = rM0MPI[i];
				p0[i + GridMPI] = pM0MPI[i];
			}
			
		}
		MPI_Bcast(r0, (NumOfGrid*NumOfGrid), MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(p0, (NumOfGrid*NumOfGrid), MPI_DOUBLE, 0, MPI_COMM_WORLD);


		#pragma omp parallel for private(i) num_threads(4)
		for (i = 0; i < GridMPI; i++) {
			r0MPI[i] = r0[i + MPICRS_ColStart];
			rM0MPI[i] = r0[i + MPICRS_RowStart];
			p0MPI[i] = p0[i + MPICRS_ColStart];
			pM0MPI[i] = p0[i + MPICRS_RowStart];
		}	

		if (SumOfN2Err < StopEps) break;
	}
	t1 = MPI_Wtime();
	
	
	FILE *myfile;
	double WritingTime1, WritingTime2, test1, test2, test3;
	// Post Reassembly
	WritingTime1 = MPI_Wtime();
	if (nrank != 0)
		MPI_Send(&x[0], GridMPI, MPI_DOUBLE, 0, 99, MPI_COMM_WORLD);
	else {
		fopen_s(&myfile, "test1file", "wb");
		fwrite(x, sizeof(double), GridMPI, myfile);
		for (i = 1; i<nproc; i++) {
			MPI_Recv(&x[0], GridMPI, MPI_DOUBLE, i, 99, MPI_COMM_WORLD, &status);
			fwrite(x, sizeof(double), GridMPI, myfile);		}
		fclose(myfile);
	}
	if (nrank == 0){
		WritingTime2 = MPI_Wtime();
		test1 = WritingTime2 - WritingTime1;}
	
	// Single task
	WritingTime1 = MPI_Wtime();
	char filename[128];
	sprintf_s(filename, 128, "test2file.%d", nrank);
	fopen_s(&myfile, filename, "wb");
	fwrite(x, sizeof(double), GridMPI, myfile);
	fclose(myfile);
	if (nrank == 0) {
		WritingTime2 = MPI_Wtime();
		test2 = WritingTime2 - WritingTime1;}

	// MPI I/O
	WritingTime1 = MPI_Wtime();
	sprintf_s(filename, 128, "test3file.%d", nrank);
	MPI_File_open(MPI_COMM_SELF, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &myfile);
	MPI_File_write(myfile, x, GridMPI, MPI_DOUBLE, MPI_STATUS_IGNORE);
	MPI_File_close(&myfile);
	if (nrank == 0) {
		WritingTime2 = MPI_Wtime();
		test3 = WritingTime2 - WritingTime1;}
	MPI_Finalize();

	// Code for debug 
	if (nrank == 0 && DebugFlag == 1) {
	double **Numerical_SolM;
	FILE *SolutionTXT;
	fopen_s(&SolutionTXT, "[MATRIX]_Solution.txt", "w");
	Numerical_SolM = (double**)malloc(sizeof(double)*NumOfGrid);
	for (i = 0; i < NumOfGrid; i++)
	Numerical_SolM[i] = (double*)malloc(sizeof(double)*NumOfGrid);
	int ik = 0, jk = 0;
	for (i = 0; i < (NumOfGrid*NumOfGrid); i++) {
	Numerical_SolM[ik][jk] = x[i];
	fprintf(SolutionTXT, "%0.20f \t", Numerical_SolM[ik][jk]);
	if ((jk - (NumOfGrid - 1)) == 0) fprintf(SolutionTXT, "\n");
	if ((i + 1) % NumOfGrid == 0) {
	ik = ik + 1; jk = 0; }
	else jk = jk + 1;}
	fclose(SolutionTXT);
	}


	if (nrank == 0) {
	double *ASol;
	double SumDiff, Error;

	ASol = (double*)malloc(sizeof(double)*NumOfGrid*NumOfGrid);

	printf("\n : Homework #1 - Parallel Scientific Computing / Junyoung Jang (2016323246)\n");
	printf("    1. Domain = [0,1]x[0,1]\n    2. # Grid = %3d\n    3. Relative Error Stop Condition = %e\n       MaxIter = %d \n\n\n", NumOfGrid, StopEps, MaxIter);
	printf("  Pure Poisson Equation DiriBD+NeumBD : [Conjugate Gradient Method]\n \t: Iteration = %6d, Time = %f, 2-Norm Error = %e, ", loop, (t1-t0), SumOfN2Err);

	ASol = AnalyticalSol_for_PurePoisson(NumOfGrid);
	
	SumDiff = 0;
	for (i = 0; i < NumOfGrid * NumOfGrid; i++)
	SumDiff += ((x[i] - ASol[i])*(x[i] - ASol[i]));

	Error = sqrt(SumDiff) / (NumOfGrid*NumOfGrid);
	printf("||e||= %e\n\n", Error);
	}


	// Post Reassembly
	if (nrank == 0) {
		printf("  [I/O] Post Reassembly : %10.5f\n", test1);
		printf("  [I/O] Single Task     : %10.5f\n", test2);
		printf("  [I/O] MPI I/O         : %10.5f\n", test3);
	}


	




	return 0;
}
////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////     Function	(END)	////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////       M A I N		////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char * argv[])
{
	int i;											// Iteration index.
	int NumOfGrid = DefineNumOfGrid;				// Number of Grid.
	double Initial_Value = 0;						// Value of Initial Matrix.
	double *x;										// Initial Value.
	double *Numerical_Sol;							// Storage of Numerical Solution.
	double *ForceTerm;								// Storage of RHS : Force Term f (Ax=f).

													// Initialize array using malloc fucntion
	x = (double*)malloc(sizeof(double)*NumOfGrid*NumOfGrid);
	ForceTerm = (double*)malloc(sizeof(double)*NumOfGrid*NumOfGrid);

	ForceTerm = ForceTerm_for_PurePoisson(NumOfGrid);

	for (i = 0; i < NumOfGrid*NumOfGrid; i++)
		x[i] = 0;

	CG_for_PurePoisson(&argc, &argv, x, ForceTerm, NumOfGrid);
}

