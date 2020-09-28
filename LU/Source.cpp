#include "LU.h"

int	main()
{
	double* A = NULL;
	double* A_copy = NULL;
	double* A_copy2 = NULL;
	double* L = NULL;
	double* U = NULL;
	double* C = NULL;
	double* b = NULL;
	double* x = NULL;

	b = new double[n];
	x = new double[n];

	double t1, t2, tms0, tms1, tms2;

	make_matrix(A);
	make_matrix(A_copy);
	make_matrix(A_copy2);
	make_matrix(L);
	make_matrix(U);
	make_matrix(C);
	init_zero(L, n * n);
	init_zero(U, n * n);
	init_zero(C, n * n);
	init_random_matrix(A);
	copy_part_matrix(A, A_copy, 0, 2, n);
	copy_part_matrix(A, A_copy2, 0, 2, n);
	//алгоритм из вики
/*	LU_decomposition(A, L, U);
	std::cout << "\nU-matrix algorithm with new matrix\n";
	print_matrix(U);
	std::cout << "\nL-matrix algorithm with new matrix\n";
	print_matrix(L);*/
	
	//решатель с помощью LU
/*	for (int i = 0; i < n; i++) b[i] = 1;
	LU_solve(L, U, b, x);
	std::cout << '\n';
	for (int i = 0; i < n; i++) std::cout << x[i] << '\t';*/

	//алгоритм из Демеля
	std::cout << "\nA-matrix algorithm with L and U instead A\n";
	t1 = omp_get_wtime();
	LU_decomposition(A);
	t2 = omp_get_wtime();
	tms0 = (double)(t2 - t1);
	std::cout << "Time LU-block = " << tms0 << " sec" << std::endl;
	std::cout << "A[n/2][n/2] = " << A[n / 2 * n + n / 2] << std::endl;
//	print_matrix(A);

	//блочное LU разложение
	std::cout << "\nLU-decomposition using blocks\n";
	t1 = omp_get_wtime();
	LU_decomposition_blocks(A_copy);
	t2 = omp_get_wtime();
	tms1 = (double)(t2 - t1);
	std::cout << "Time LU-block = " << tms1 << " sec" << std::endl;
	std::cout << "A[n/2][n/2] = " << A_copy[n / 2 * n + n / 2] << std::endl;
//	print_matrix(A_copy);

	std::cout << "Times with usual U = " << tms0 / tms1 << std::endl;

	//блочное LU разложение (транспонированная U)
	std::cout << "\nLU-decomposition using blocks with U transpose\n";
	t1 = omp_get_wtime();
	LU_transpose_decomposition_blocks(A_copy2);
	t2 = omp_get_wtime();
	tms2 = (double)(t2 - t1);
	std::cout << "Time LU-block = " << tms2 << " sec" << std::endl;
	std::cout << "A[n/2][n/2] = " << A_copy2[n / 2 * n + n / 2] << std::endl;
	std::cout << "Times with U transpose = " << tms0 / tms1 << std::endl;
//		print_matrix(A_copy2);

	//очистка памяти
	clear_matrix(A);
	clear_matrix(L);
	clear_matrix(U);
	clear_matrix(C);
	return (0);
}