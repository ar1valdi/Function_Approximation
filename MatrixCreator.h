#pragma once
#include "Matrix.h"

class MatrixCreator {
private:	
	static void init_edge_cases_diag_5(Matrix& m, int, double, double, double);
public:
	static Matrix diagonal_5(int N, double a1, double a2, double a3);
	static Matrix diagonal_3(int N, double a1, double a2);
	static Matrix diagonal_1(int N, double a1);
	static Matrix vector_sin_n_f_1(int N);
	static Matrix vector_n(int N);
};