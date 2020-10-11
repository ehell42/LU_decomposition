#pragma once
#include <iostream>
#include <fstream>
#include "omp.h"

const int n = 1024;	//размерность матрицы
const int b = 64;	//размерность блока

void	make_matrix(double* &A);
void	init_matrix(double*& A);
void	init_random_matrix(double*& A);
void	init_zero(double*& A, int num);
void	LU_decomposition(double*& A, double*& L, double*& U);
void	LU_decomposition(double*& A);
void	LU_decomposition(double*& A, int m);
void	LU_decomposition_blocks(double*& A);
void	LU_solve(double* L, double* U, double* b, double*& x);
void	print_matrix(double* A);
void	mult(double* &C, double *L, double *U);
void	clear_matrix(double* &A);
void	copy_part_matrix(double* A, double* &B, int d, int f, int m);
void	LU_transpose_decomposition_blocks(double*& A);
void	LU_transpose_decomposition_blocks_parallel(double*& A);