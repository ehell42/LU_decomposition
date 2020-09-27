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

	file.open("matrix.txt");
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