#include "Matrix.h"
#include <cmath>
#include <iomanip>

Matrix::Matrix(){
	n = m = 0;
}
Matrix::Matrix(double n) {
	this->n = m = n;
	mat.resize(n, vector<double>(n, 0));
}
Matrix::Matrix(double n, double m) {
	this->n = n;
	this->m = m;
	mat.resize(n, vector<double>(m, 0));
}
Matrix::Matrix(double n, double m, double x) {
	this->n = n;
	this->m = m;
	mat.resize(n, vector<double>(m, x));
}
Matrix::Matrix(vector<vector<double>> v) {
	n = v.size();
	m = v[0].size();
	mat = v;
}
void Matrix::print() const {
	cout << std::fixed << std::setprecision(3);
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
            cout << std::setw(6) << mat[j][i] << "  ";
		}
		putchar('\n');
	}
	cout << std::scientific;
}
pair<int, int> Matrix::getSize() const {
	return { n, m };
}
double Matrix::get_norm() const {
	double norm = 0;

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			norm += mat[i][j] * mat[i][j];
		}
	}
	return sqrt(norm);
}
const Matrix& Matrix::operator=(const Matrix& m1) {
	mat = m1.mat;
	n = m1.n;
	m = m1.m;
	return *this;
}