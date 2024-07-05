#pragma once
#include "Matrix.h"
#include <ctime>
#include <list>

struct ExReturnData {
	Matrix m;
	std::list<double> errors;
	int iters;
	double time;
};