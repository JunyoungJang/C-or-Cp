// By Junyoung Jang - Ph.D Student of Yonsei CSE (Updated Date : 2016. 09. 26)
#include "pde.h"

/////////////////////////////////////////////////////////////////////////////////////////////////////
/* Iterative Method for Pure Poisson Equation.
Pure Poisson Equation : L(u)=f where L is Laplace Operator
1. Jacobi Method					- Jacobi_for_PurePoisson		(Updated Date : 2016. 09. 24.)
2. Gauss-Seidel Method with SOR		- GaussSOR_for_PurePoisson		(Updated Date : 2016. 09. 24.)
3. Conjugate Gradient Method		- CG_for_PurePoisson			(Updated Date : 2016. 09. 24.)*/
/////////////////////////////////////////////////////////////////////////////////////////////////////
double(*AnalyticalSol_for_PurePoisson(double ASol[MatrixSize][MatrixSize], int NumOfGrid))[MatrixSize] {
	int i, j;											// Iteration Index.
	double dxdygrid = (NumOfGrid - 1);
	FILE *SolutionTXT;

	fopen_s(&SolutionTXT, "[MATRIX]_AnalyticalSolution.txt", "w");
	for (i = 0;i < NumOfGrid;i++) {
		for (j = 0;j < NumOfGrid;j++) {
			ASol[i][j] = -0.5*(sin(PI*i/dxdygrid)*cos(PI*j/dxdygrid) / (PI*PI));
			fprintf(SolutionTXT, "%0.20f \t", ASol[i][j]);
			if ((j - (NumOfGrid - 1)) == 0) fprintf(SolutionTXT, "\n");}}
	fclose(SolutionTXT);
	return ASol;
}

double(*ForceTerm_for_PurePoisson(double Fvalue[MatrixSize][MatrixSize], int NumOfGrid))[MatrixSize] {
	int i, j;											// Iteration Index.
	double dxdygrid = (NumOfGrid - 1);

	for (i = 0;i < NumOfGrid;i++) {
		for (j = 0;j < NumOfGrid;j++) {
			Fvalue[i][j] = (sin(PI*i/ dxdygrid)*cos(PI*j/ dxdygrid)) / ((dxdygrid+1)*(dxdygrid+1));}}
	return Fvalue;
}

double(*Jacobi_for_PurePoisson(double x0[MatrixSize][MatrixSize], double Fvalue[MatrixSize][MatrixSize], int NumOfGrid))[MatrixSize] {
	int i, j;											// Iteration Index.
	int k = 0;											// #Iteration.
	double dxdygrid = (NumOfGrid - 1);
	double N2Err = 0;									// (Sub) Stop Condition.
	double SumOfN2Err = 1;								// Initial Stop Condition for (While).
	double new_x[MatrixSize][MatrixSize] = { 0 };		// Storage of Updated Matrix
	clock_t start, end;
	FILE *RelErrorTXT, *SolutionTXT;

	start = clock();
	fopen_s(&RelErrorTXT, "[RelError]_Jacobi.txt", "w");
	while (SumOfN2Err > StopEps) {
		// Core of Jacobi Method.
		for (j = 1;j < (NumOfGrid-1);j++) {
			for (i = 1; i < (NumOfGrid-1);i++) {
				new_x[i][0] = ((x0[i - 1][0] + x0[i + 1][0] + 2 * x0[i][1]) - Fvalue[i][0]) * 0.25;
				new_x[i][j] = ((x0[i - 1][j] + x0[i + 1][j] + x0[i][j - 1] + x0[i][j + 1]) - Fvalue[i][j]) * 0.25;
				new_x[i][NumOfGrid - 1] = ((x0[i - 1][NumOfGrid - 1] + x0[i + 1][NumOfGrid - 1] 
													+ 2 * x0[i][NumOfGrid - 2]) - Fvalue[i][NumOfGrid - 1]) * 0.25;}}

		for (j = 0;j < (NumOfGrid - 1);j++) {
			new_x[0][j] = 0;new_x[NumOfGrid - 1][j] = 0;}

		// 2-Norm Error (Stop Condition) of Iterative Method.
		SumOfN2Err = 0;
		for (j = 0;j < (NumOfGrid);j++) {
			for (i = 0; i < (NumOfGrid);i++) {
				N2Err = (new_x[i][j] - x0[i][j])*(new_x[i][j] - x0[i][j]);
				SumOfN2Err += N2Err;}}
		SumOfN2Err = sqrt(SumOfN2Err);
		fprintf(RelErrorTXT, "%6d \t %0.20f \n", k, SumOfN2Err);
		
		// Update Recent Matrix (x0) to Previous Matrix (new_x).
		for (j = 0;j < (NumOfGrid);j++) {
			for (i = 0; i < (NumOfGrid);i++) {
				x0[i][j] = new_x[i][j];}}
		if (k == MaxIter) break;
		k = k + 1;}

	fclose(RelErrorTXT);
	end = clock();
	printf("  Pure Poisson Equation : [Jacobi Method]\n \t: Iteration = %6d, Time = %f, 2-Norm Error = %0.20f, ", k, (double)(end - start) / CLOCKS_PER_SEC, SumOfN2Err);

	fopen_s(&SolutionTXT, "[MATRIX]_Jacobi.txt", "w");
	for (i = 0;i < NumOfGrid;i++) {
		for (j = 0;j < NumOfGrid;j++) {
			fprintf(SolutionTXT, "%0.20f \t", x0[i][j]);
			if ((j - (NumOfGrid - 1)) == 0) fprintf(SolutionTXT, "\n");}}
	fclose(SolutionTXT);

	return x0;
}

double(*GaussSOR_for_PurePoisson(double x0[MatrixSize][MatrixSize], double Fvalue[MatrixSize][MatrixSize], int NumOfGrid))[MatrixSize] {
	int i, j;												// Iteration Index.
	int k = 0;												// #Iteration.
	double dxdygrid = (NumOfGrid - 1);					
	double N2Err = 0;										// (Sub) Stop Condition.
	double SumOfN2Err = 1;									// Initial Stop Condition for (While).
	double new_x[MatrixSize][MatrixSize] = { 0 };			// Storage of Updated Matrix 
	double OptWeight = 2 / (1 + PI / NumOfGrid);			// Optimal Weight(omega) of SOR.
	clock_t start, end;
	FILE *RelErrorTXT;
	FILE *SolutionTXT;
	start = clock();
	fopen_s(&RelErrorTXT, "[RelError]_GaussWithSOR.txt", "w");

	while (SumOfN2Err > StopEps) {
		// Update Previous Matrix (x0) to Recent Matrix (new_x).
		for (j = 0;j < (NumOfGrid);j++) {
			for (i = 0; i < (NumOfGrid);i++) {
				new_x[i][j] = x0[i][j];}}

		// Core of Gauss with SOR Method.
		for (j = 1;j < (NumOfGrid - 1);j++) {
			for (i = 1; i < (NumOfGrid - 1);i++) {
				x0[i][0] = (1 - OptWeight)*x0[i][0] +
					OptWeight*(0.25*((x0[i - 1][0] + x0[i + 1][0] + 2*x0[i][1]) - Fvalue[i][0]));
				x0[i][j] = (1 - OptWeight)*x0[i][j] +
					OptWeight*(0.25*((x0[i - 1][j] + x0[i + 1][j] + x0[i][j - 1] + x0[i][j + 1]) - Fvalue[i][j]));
				x0[i][NumOfGrid - 1] = (1 - OptWeight)*x0[i][NumOfGrid - 1] +
					OptWeight*0.25*((x0[i - 1][NumOfGrid - 1] + x0[i + 1][NumOfGrid - 1] 
												+ 2*x0[i][NumOfGrid - 2] - Fvalue[i][NumOfGrid - 1]));}}
	
		for (j = 0;j < (NumOfGrid);j++) {
			x0[0][j] = 0;
			x0[NumOfGrid - 1][j] = 0;}

		// 2-Norm Error (Stop Condition) of Iterative Method.
		SumOfN2Err = 0;
		for (j = 0;j < (NumOfGrid);j++) {
			for (i = 0; i < (NumOfGrid);i++) {
				N2Err = (new_x[i][j] - x0[i][j])*(new_x[i][j] - x0[i][j]);
				SumOfN2Err += N2Err;}}
		SumOfN2Err = sqrt(SumOfN2Err);
		fprintf(RelErrorTXT, "%6d \t %0.20f \n", k, SumOfN2Err);

		if (k == MaxIter) break;
		k = k + 1;}

	fclose(RelErrorTXT);
	end = clock();
	printf("  Pure Poisson Equation : [Gauss(SOR) Method]\n \t: Iteration = %6d, Time = %f, 2-Norm Error = %0.20f, ", k, (double)(end - start) / CLOCKS_PER_SEC, SumOfN2Err);

	fopen_s(&SolutionTXT, "[MATRIX]_GaussWithSOR.txt", "w");
	for (i = 0;i < NumOfGrid;i++) {
		for (j = 0;j < NumOfGrid;j++) {
			fprintf(SolutionTXT, "%0.20f \t", x0[i][j]);
			if ((j - (NumOfGrid - 1)) == 0) fprintf(SolutionTXT, "\n");}}
	fclose(SolutionTXT);

	return x0;
}

double(*CG_for_PurePoisson(double x0[MatrixSize][MatrixSize], double Fvaluef[MatrixSize][MatrixSize], int NumOfGrid))[MatrixSize] {
	int i, j, ik, jk;											// Iteration Index.
	int ITER = 0;											// #Iteration.
	double dxdygrid = (NumOfGrid - 1);
	double N2Err = 0;										// (Sub) Stop Condition.
	double SumOfN2Err = 2;									// Initial Stop Condition for (While).
	double x[MatrixSize*MatrixSize] = { 0 };			// Storage of Updated Matrix 
	double new_x[MatrixSize*MatrixSize] = { 0 };			// Storage of Updated Matrix 
	double r0[MatrixSize*MatrixSize] = { 0 };				// Storage of OLD Residual Matrix.
	double p0[MatrixSize*MatrixSize] = { 0 };				// Storage of OLD Direction Matrix. 
	double FvalueVector[MatrixSize*MatrixSize] = { 0 };
	double PoissonMultipleVector[MatrixSize*MatrixSize] = { 0 };			// Storage of Updated Matrix 
	double Alpha, Beta;										// CG Parameter.
	double UpperValue, LowerValue;
	clock_t start, end;
	FILE *RelErrorTXT;
	FILE *SolutionTXT;
	start = clock();
	fopen_s(&RelErrorTXT, "[RelError]_Conjugate_Gradient.txt", "w");

	// ForceTerm and Initial Conditon Vectorization. 
	for (i = 0;i < NumOfGrid;i++) {
		for (j = 0;j < NumOfGrid;j++) {
			FvalueVector[ITER] = Fvaluef[i][j];
			x[ITER] = x0[i][j];
			ITER = ITER + 1;}}
	ITER = 0;

	// Construct Pure Poisson Matrix (Matrix-free) 
	PoissonMultipleVector[0] = -4 * x[0] + x[1] + x[NumOfGrid];
	for (i = 1;i < (NumOfGrid*NumOfGrid - 1);i++) {
		PoissonMultipleVector[i] = x[i - 1] - 4 * x[i] + x[i + 1] + x[i + NumOfGrid];
		if ((i + 1) % (NumOfGrid) == 0) PoissonMultipleVector[i] = PoissonMultipleVector[i] - x[i + 1];
		if (i % (NumOfGrid) == 0) PoissonMultipleVector[i] = PoissonMultipleVector[i] - x[i - 1];
		if (i > (NumOfGrid - 1)) PoissonMultipleVector[i] = PoissonMultipleVector[i] + x[i - (NumOfGrid)];}
	PoissonMultipleVector[NumOfGrid*NumOfGrid - 1] = x[NumOfGrid*NumOfGrid - (NumOfGrid + 1)] - 4 * x[NumOfGrid*NumOfGrid - 1] + x[NumOfGrid*NumOfGrid - 2];

	// Set Initial Value.
	for (i = 0;i < (NumOfGrid*NumOfGrid);i++) {
		r0[i] = -PoissonMultipleVector[i] + FvalueVector[i];
		p0[i] = r0[i];
	}

	while (SumOfN2Err > 0.001) {

		for (i = 0; i < (NumOfGrid*NumOfGrid); i++) {
			x[i] = new_x[i];
		}


		// Calculate Conjugate Parameter Alpha.
		UpperValue = 0;
		for (i = 0; i < (NumOfGrid*NumOfGrid); i++) {
			UpperValue += p0[i] * r0[i];
		}

		for (i = 0;i < (NumOfGrid*NumOfGrid);i++) {						// Construct Pure Poisson Matrix (Matrix-free) 
			PoissonMultipleVector[i] = 0;
		}
		PoissonMultipleVector[0] = -4 * p0[0] + p0[1] + p0[NumOfGrid];
		for (i = 1;i < (NumOfGrid*NumOfGrid - 1);i++) {
			PoissonMultipleVector[i] = p0[i - 1] - 4 * p0[i] + x[i + 1] + x[i + NumOfGrid];
			if ((i + 1) % (NumOfGrid) == 0) PoissonMultipleVector[i] = PoissonMultipleVector[i] - p0[i + 1];
			if (i % (NumOfGrid) == 0) PoissonMultipleVector[i] = PoissonMultipleVector[i] - p0[i - 1];
			if (i > (NumOfGrid - 1)) PoissonMultipleVector[i] = PoissonMultipleVector[i] + p0[i - (NumOfGrid)];
		}
		PoissonMultipleVector[NumOfGrid*NumOfGrid - 1] =
			x[NumOfGrid*NumOfGrid - (NumOfGrid + 1)] - 4 * x[NumOfGrid*NumOfGrid - 1] + x[NumOfGrid*NumOfGrid - 2];

		LowerValue = 0;
		for (i = 0; i < (NumOfGrid*NumOfGrid); i++) {
			LowerValue += p0[i] * PoissonMultipleVector[i];
		}
		Alpha = UpperValue / LowerValue;


		// Calculate Optimal Vector
		for (i = 0;i < (NumOfGrid*NumOfGrid);i++) {
			new_x[i] = x[i] + Alpha * p0[i];
			r0[i] = r0[i] - Alpha * PoissonMultipleVector[i];
		}

		// Calculate Conjugate Parameter Beta.
		UpperValue = 0;
		for (i = 0;i < (NumOfGrid*NumOfGrid);i++) {
			UpperValue += r0[i] * PoissonMultipleVector[i];
		}
		Beta = -UpperValue / LowerValue;

		// Calculate Conjugate Parameter Beta.
		for (i = 0;i < (NumOfGrid*NumOfGrid);i++) {
			p0[i] = r0[i] - Beta * p0[i];
		}
		SumOfN2Err = 0;
		for (i = 0;i < (NumOfGrid*NumOfGrid);i++) {
			N2Err = (new_x[i] - x[i])*(new_x[i] - x[i]);
			SumOfN2Err += N2Err;
		}
		SumOfN2Err = sqrt(SumOfN2Err);
		ITER = ITER + 1;
	}
	fclose(RelErrorTXT);
	end = clock();
	printf("  Pure Poisson Equation : [Conjugate Gradient Method]\n \t: Iteration = %6d, Time = %f, 2-Norm Error = %0.20f, ", ITER, (double)(end - start) / CLOCKS_PER_SEC, SumOfN2Err);

	ik = 0;jk = 0;
	for (i = 0;i < (NumOfGrid*NumOfGrid);i++) {
		x0[ik][jk] = x[i];
		if ((i + 1) % NumOfGrid == 0){
			ik = ik + 1;jk = 0;
		}
		else jk = jk + 1;
	}

	fopen_s(&SolutionTXT, "[MATRIX]_Conjugate_Gradient.txt", "w");
	for (i = 0;i < NumOfGrid;i++) {
		for (j = 0;j < NumOfGrid;j++) {
		fprintf(SolutionTXT, "%0.20f \t", x0[i][j]);
			if ((j - (NumOfGrid - 1)) == 0) fprintf(SolutionTXT, "\n");}}
	fclose(SolutionTXT);
	return x0;
}

	
	
	


	

