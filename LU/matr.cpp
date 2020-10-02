#include "LU.h"

void	make_matrix(double*& A)
{
	A = new double[n * n];
}
void	clear_matrix(double*& A)
{
	delete[] A;
}
void	init_matrix(double*& A)
{
	std::ifstream file;

	file.open("rand_matr2048.txt");
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			file >> A[i * n + j];
	file.close();
}
void	init_random_matrix(double*& A)
{
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			A[i * n + j] = rand() % n + 1.0;
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
void	init_zero(double*& A, int num)
{
	for (int i = 0; i < num; i++)
		A[i] = 0;
}
//works right
void	copy_part_matrix(double* A, double*& B, int d, int f, int m)
{
	if (f == 0)
	{
		for (int i = d - m; i < d; i++)
			for (int j = d; j < n; j++)
				B[(i - (d - m)) * (n - d) + (j - d)] = A[i * n + j];
	}
	else if (f == 1)
	{
		for (int i = d; i < n; i++)
			for (int j = d - m; j < d; j++)
				B[(i - d) * m + (j - (d - m))] = A[i * n + j];
	}
	else if (f == 2)
	{
		for (int i = d; i < d + m; i++)
			for (int j = d; j < d + m; j++)
				B[(i - d) * m + (j - d)] = A[i * n + j];
	}
	else
	{
		for (int i = d - m; i < d; i++)
			for (int j = d; j < n; j++)
				B[(j - d) * m + i - (d - m)] = A[i * n + j];
	}
}