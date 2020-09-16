#include <iostream>

const int n = 5;

void	make_matrix(double* &A);
void	init_matrix(double*& A);
void	init_zero(double*& A);
void	LU_decomposition(double*& A, double*& L, double*& U);
void	LU_decomposition(double*& A);
void	print_matrix(double* A);
void	mult(double* &C, double *L, double *U);

int	main()
{
	double* A = NULL;
	double* L = NULL;
	double* U = NULL;
	double* C = NULL;

	make_matrix(A);
	make_matrix(L);
	make_matrix(U);
	make_matrix(C);
	init_matrix(A);
	init_zero(L);
	init_zero(U);
	init_zero(C);
	std::cout << "A-matrix:\n";
	print_matrix(A);

	//алгоритм из вики
	LU_decomposition(A, L, U);
	std::cout << "\nU-matrix algorithm with new matrix\n";
	print_matrix(U);
	std::cout << "\nL-matrix algorithm with new matrix\n";
	print_matrix(L);

	//алгоритм из Демеля
	std::cout << "\nNew A-matrix algorithm with L and U instead A\n";
	LU_decomposition(A);
	print_matrix(A);
	return (0);
}

void	make_matrix(double*& A)
{
	A = new double[n * n];
}

void	init_matrix(double*& A)
{
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			A[i * n + j] = sin(i + 0.5) * cos(j+0.5);
}

void	init_zero(double*& A)
{
	for (int i = 0; i < n * n; i++)
		A[i] = 0;
}

void	LU_decomposition(double*& A, double*& L, double*& U)
{
	for (int j = 0; j < n; j++)
	{
		U[j] = A[j];
		L[j * n] = A[j * n] / U[0];
	}//print_matrix(L); print_matrix(U);
	for (int i = 1; i < n; i++)
	{
		for (int j = i; j < n; j++)
		{
			U[i * n + j] = A[i * n + j];
			for (int k = 0; k < i - 1; k++)
				U[i * n + j] -= L[i * n + k] * U[k * n + j];

			L[j * n + i] = 1. / U[i * n + i] * A[j * n + i];// std::cout << U[i * n + i] << '\n';
			for (int k = 0; k < i - 1; k++)
				L[j * n + i] -= 1. / U[i * n + i] * L[j * n + k] * U[k * n + i];
		}
	}
}

void	LU_decomposition(double*& A)
{
	for (int i = 0; i < n - 1; i++)
	{
		for (int j = i + 1; j < n; j++)
			A[j * n + i] = A[j * n + i] / A[i * n + i];
		for (int j = i + 1; i < n; i++)
			for (int k = i + 1; k < n; k++)
				A[j * n + k] = A[j * n + k] - A[j * n + i] * A[i * n + k];
	}
}

void	print_matrix(double* A)
{
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++)
			std::cout << A[i * n + j] << "\t";
		std::cout << '\n';
	}
}

void	mult(double*& C, double* L, double* U)
{
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			for (int k = 0; k < n; k++)
				C[i * n + j] += L[i * n + k] * U[k * n + j];
}