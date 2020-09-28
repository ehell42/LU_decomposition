#include "LU.h"

//from wiki (works)
void	LU_decomposition(double*& A, double*& L, double*& U)
{
	for (int j = 0; j < n; j++)
	{
		U[j] = A[j];
		L[j * n] = A[j * n] / U[0];
	}
	for (int i = 1; i < n; i++)
	{
		for (int j = i; j < n; j++)
		{
			U[i * n + j] = A[i * n + j];
			for (int k = 0; k < i; k++)
				U[i * n + j] -= L[i * n + k] * U[k * n + j];
			L[j * n + i] = 1. / U[i * n + i] * A[j * n + i];
			for (int k = 0; k < i; k++)
				L[j * n + i] -= 1. / U[i * n + i] * L[j * n + k] * U[k * n + i];
		}
	}
}
//from Demmel (works)
void	LU_decomposition(double*& A)
{
	for (int i = 0; i < n - 1; i++)
	{
		for (int j = i + 1; j < n; j++)
			A[j * n + i] = A[j * n + i] / A[i * n + i];
		for (int j = i + 1; j < n; j++)
			for (int k = i + 1; k < n; k++)
				A[j * n + k] = A[j * n + k] - A[j * n + i] * A[i * n + k];
	}
}
void	LU_decomposition(double*& A, int m)
{
	for (int i = 0; i < m - 1; i++)
	{
		for (int j = i + 1; j < m; j++)
			A[j * m + i] = A[j * m + i] / A[i * m + i];
		for (int j = i + 1; j < m; j++)
			for (int k = i + 1; k < m; k++)
				A[j * m + k] = A[j * m + k] - A[j * m + i] * A[i * m + k];
	}
}
//blocks LU-decomposition from Demmel (works!)
void	LU_decomposition_blocks(double*& A)
{
	double* A12_tmp = new double[(n - b) * b];
	double* A21_tmp = new double[(n - b) * b];
	double* A11_tmp = new double[b * b];

	init_zero(A12_tmp, (n - b) * b);
	init_zero(A21_tmp, (n - b) * b);
	init_zero(A11_tmp, b * b);

	for (int i = 0; i < n; i += b)
	{
		//find diagonal block (works)
		copy_part_matrix(A, A11_tmp, i, 2, b);//right
		LU_decomposition(A11_tmp, b);//right
		//copy side matrix
		copy_part_matrix(A, A12_tmp, i + b, 0, b);//right
		copy_part_matrix(A, A21_tmp, i + b, 1, b);//right
		//find U_12 (works)
		for (int j = 0; j < n - (i + b); j++)
			for (int k = 1; k < b; k++)
				for (int l = 0; l < k; l++)
					A12_tmp[k * (n - (i + b)) + j] -= A12_tmp[l * (n - (i + b)) + j] * A11_tmp[k * b + l];
		//find L_21 (works)
		for (int j = 0; j < n - (i + b); j++)
		{
			A21_tmp[j * b] /= A11_tmp[0];
			for (int k = 1; k < b; k++)
			{
				for (int l = 0; l < k; l++)
					A21_tmp[j * b + k] -= A21_tmp[j * b + l] * A11_tmp[l * b + k];
				A21_tmp[j * b + k] /= A11_tmp[k * b + k];
			}
		}
		//change matrix A22 (works)
		for (int j = i + b; j < n; j++)
			for (int l = i + b; l < n; l++)
				for (int k = 0; k < b; k++)
					A[j * n + l] -= A21_tmp[(j - (i + b)) * b + k] * A12_tmp[k * (n - (b + i)) + (l - (i + b))];
		//copy matrix A11 in A (works)
		for (int j = i; j < i + b; j++)
			for (int l = i; l < i + b; l++)
				A[j * n + l] = A11_tmp[(j - i) * b + l - i];
		//copy matrix A12 in A (works)
		for (int j = i; j < b + i; j++)
			for (int l = i + b; l < n; l++)
				A[j * n + l] = A12_tmp[(j - i) * (n - (i + b)) + l - (i + b)];
		//copy matrix A21 in A (works)
		for (int j = i + b; j < n; j++)
			for (int l = i; l < i + b; l++)
				A[j * n + l] = A21_tmp[(j - (i + b)) * b + l - i];
	}
	clear_matrix(A11_tmp);
	clear_matrix(A21_tmp);
	clear_matrix(A12_tmp);
}
//blocks LU-decomposition from Demmel with transpose U
void	LU_transpose_decomposition_blocks(double*& A)
{
	double* A12_tmp = new double[(n - b) * b];
	double* A21_tmp = new double[(n - b) * b];
	double* A11_tmp = new double[b * b];

	init_zero(A12_tmp, (n - b) * b);
	init_zero(A21_tmp, (n - b) * b);
	init_zero(A11_tmp, b * b);

	for (int i = 0; i < n; i += b)
	{
		//find diagonal block (works)
		copy_part_matrix(A, A11_tmp, i, 2, b);//right
		LU_decomposition(A11_tmp, b);//right
		//copy side matrix
		copy_part_matrix(A, A12_tmp, i + b, 3, b);//right
		copy_part_matrix(A, A21_tmp, i + b, 1, b);//right
		//find U_12 (works)
		for (int j = 0; j < n - (i + b); j++)
			for (int k = 1; k < b; k++)
				for (int l = 0; l < k; l++)
					A12_tmp[j * b + k] -= A12_tmp[j * b + l] * A11_tmp[k * b + l];
		//find L_21 (works)
		for (int j = 0; j < n - (i + b); j++)
		{
			A21_tmp[j * b] /= A11_tmp[0];
			for (int k = 1; k < b; k++)
			{
				for (int l = 0; l < k; l++)
					A21_tmp[j * b + k] -= A21_tmp[j * b + l] * A11_tmp[l * b + k];
				A21_tmp[j * b + k] /= A11_tmp[k * b + k];
			}
		}
		//change matrix A22 (works)
		for (int j = i + b; j < n; j++)
			for (int l = i + b; l < n; l++)
				for (int k = 0; k < b; k++)
					A[j * n + l] -= A21_tmp[(j - (i + b)) * b + k] * A12_tmp[(l - (i + b)) * b + k];
		//copy matrix A11 in A (works)
		for (int j = i; j < i + b; j++)
			for (int l = i; l < i + b; l++)
				A[j * n + l] = A11_tmp[(j - i) * b + l - i];
		//copy matrix A12 in A (works)
		for (int j = i; j < b + i; j++)
			for (int l = i + b; l < n; l++)
				A[j * n + l] = A12_tmp[(l - (i + b)) * b + j - i];
		//copy matrix A21 in A (works)
		for (int j = i + b; j < n; j++)
			for (int l = i; l < i + b; l++)
				A[j * n + l] = A21_tmp[(j - (i + b)) * b + l - i];
	}
	clear_matrix(A11_tmp);
	clear_matrix(A21_tmp);
	clear_matrix(A12_tmp);
}

//solver (works)
void	LU_solve(double* L, double* U, double* b, double*& x)
{
	double* y = new double[n];

	//прямой ход Гаусса для нахождения y (Ly=b)
	y[0] = b[0] / L[0];
	for (int i = 1; i < n; i++)
	{
		y[i] = b[i] / L[i * n + i];
		for (int k = 0; k < i; k++)
			y[i] -= y[k] * L[i * n + k];
	}
	//обратный ход Гаусса для нахождения x (Ux = y)
	x[n - 1] = y[n - 1] / U[(n - 1) * n + (n - 1)];
	for (int i = n - 2; i >= 0; i--)
	{
		x[i] = y[i] / U[i * n + i];
		for (int k = i + 1; k < n; k++)
			x[i] -= x[k] * U[i * n + k];
	}
	delete[] y;
}