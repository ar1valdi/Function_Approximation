#pragma once
#include <list>

//193237
namespace Consts {
	const int C = 3;
	const int D = 7;
	const int E = 2;
	const int F = 3;
	const int N = 900 + 10 * C + D;

	const double A1_1 = 5 + E;
	const double A2_1 = -1;
	const double A3_1 = -1;

	const double A1_2 = 3;
	const double A2_2 = -1;
	const double A3_2 = -1;

	const int TEST_SIZE = 50;
	const float RESIDUUM_ACCEPTANCE = 1E-9;
	const int MAX_ITER = 100;

	const int JACOBI = 1;
	const int GAUSS_SEIDEL = 2;

	const int SCATTER_SIZE = 40;

	const std::list<int> LENGHTS_TO_COMPARE = { 100, 500, 1000, 2000, 3000};
};

