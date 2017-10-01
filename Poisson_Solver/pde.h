// By Junyoung Jang - Ph.D Student of Yonsei CSE (Updated Date : 2016. 09. 26)
#include <math.h>
#include <stdio.h>
#include <time.h>

#define PI 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651328
#define StopEps 1e-5							// Stop Condition.
#define MatrixSize 100						    // Maximum #Grid.
#define MaxIter 10000							// Maximum #Iteration. 
#define DefineNumOfGrid 50

/////////////////////////////////////////////////////////////////////////////////////////////////////
/* Iterative Method for Pure Poisson Equation.
Pure Poisson Equation : L(u)=f where L is Laplace Operator
1. Jacobi Method					- Jacobi_for_PurePoisson()		(Updated Date : 2016. 09. 24.)
2. Gauss-Seidel Method with SOR		- GaussSOR_for_PurePoisson()	(Updated Date : 2016. 09. 24.)
3. Conjugate Gradient Method		- CG_for_PurePoisson()			(Updated Date : 2016. 09. 24.)*/
/////////////////////////////////////////////////////////////////////////////////////////////////////

double(*AnalyticalSol_for_PurePoisson(double ASol[MatrixSize][MatrixSize], int NumOfGrid))[MatrixSize];
double(*ForceTerm_for_PurePoisson(double Fvalue[MatrixSize][MatrixSize], int NumOfGrid))[MatrixSize];
double(*Jacobi_for_PurePoisson(double x0[MatrixSize][MatrixSize], double f[MatrixSize][MatrixSize], int NumOfGrid))[MatrixSize];
double(*GaussSOR_for_PurePoisson(double x0[MatrixSize][MatrixSize], double f[MatrixSize][MatrixSize], int NumOfGrid))[MatrixSize];
double(*CG_for_PurePoisson(double x0[MatrixSize][MatrixSize], double Fvaluef[MatrixSize][MatrixSize], int NumOfGrid))[MatrixSize];