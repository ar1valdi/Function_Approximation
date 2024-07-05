#pragma once
#include "Matrix.h"

Matrix operator+(const Matrix& m1, const Matrix& m2) {
	pair<int, int> size = m1.getSize();
	if (size != m2.getSize()) {
		throw new runtime_error("incompatible sizes");
	}
	Matrix m = m1;
	for (int i = 0; i < size.first; i++) {
		for (int j = 0; j < size.second; j++) {
			m[i][j] = m1[i][j] + m2[i][j];
		}
	}
	return m;
}

Matrix operator-(const Matrix& m1, const Matrix& m2) {
	pair<int, int> size = m1.getSize();
	if (size != m2.getSize()) {
		throw new runtime_error("incompatible sizes");
	}
	Matrix m = m1;
	for (int i = 0; i < size.first; i++) {
		for (int j = 0; j < size.second; j++) {
			m[i][j] = m1[i][j] - m2[i][j];
		}
	}
	return m;
}

Matrix operator*(const Matrix& m1, const Matrix& m2) {
	pair<int, int> size1 = m1.getSize();
	pair<int, int> size2 = m2.getSize();
	if (size1.first != size2.second) {
		throw new runtime_error("incompatible sizes");
	}
	Matrix m(size2.first, size1.second);
	for (int i = 0; i < size2.first; i++) {
		for (int j = 0; j < size1.second; j++) {
			for (int k = 0; k < size2.second; k++) {
				m[i][j] += m1[k][j] * m2[i][k];
			}
		}
	}
	return m;
}

vector<double> operator*(const vector<double>& v, double d) {
	vector<double> res = v;
	for (int i = 0; i < res.size(); i++) {
		res[i] *= d;
	}
	return res;
}

vector<double> operator-(const vector<double>& v1, const vector<double>& v2) {
	vector<double> res = v1;
	for (int i = 0; i < res.size(); i++) {
		res[i] -= v2[i];
	}
	return res;
}