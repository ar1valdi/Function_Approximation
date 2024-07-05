#pragma once
#include <iostream>
#include <vector>
using namespace std;

class Matrix {
private:
	int n, m;
	vector<vector<double>> mat;
public:
	Matrix();
	Matrix(double n);
	Matrix(double n, double m);
	Matrix(double n, double m, double x);
	Matrix(vector<vector<double>>);

	void print() const;
	vector<double>& operator[](int index) {
		return mat[index];
	}
	const vector<double>& operator[](int index) const {
		return mat[index];
	}
	const Matrix& operator*(double d) {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				mat[i][j] *= d;
			}
		}
	}
	pair<int, int> getSize() const;
	double get_norm() const;
	const Matrix& operator=(const Matrix&);
};