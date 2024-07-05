#include "MatrixCreator.h"
#include <cmath>
#include "Consts.h"

Matrix MatrixCreator::diagonal_5(int N, double a1, double a2, double a3) {
    Matrix mat(N);
    init_edge_cases_diag_5(mat, N, a1, a2, a3);
    for (int i = 2; i < N - 2; i++) {
        mat[i][i - 2] = a3;
        mat[i][i - 1] = a2;
        mat[i][i] = a1;
        mat[i][i + 1] = a2;
        mat[i][i + 2] = a3;
    }
    return mat;
}
void MatrixCreator::init_edge_cases_diag_5(Matrix& mat, int N, double a1, double a2, double a3) {
    mat[0][0] = a1;
    mat[0][1] = a2;
    mat[0][2] = a3;

    mat[1][0] = a2;
    mat[1][1] = a1;
    mat[1][2] = a2;
    mat[1][3] = a3;

    mat[N - 2][N - 4] = a3;
    mat[N - 2][N - 3] = a2;
    mat[N - 2][N - 2] = a1;
    mat[N - 2][N - 1] = a2;

    mat[N - 1][N - 3] = a3;
    mat[N - 1][N - 2] = a2;
    mat[N - 1][N - 1] = a1;
}
Matrix MatrixCreator::vector_sin_n_f_1(int N) {
    Matrix m(1, N);
    for (int i = 0; i < N; i++) {
        m[0][i] = sin((i + 1)*(Consts::F + 1));
    }
    return m;
}
Matrix MatrixCreator::vector_n(int N) {
    Matrix m(1, N);
    for (int i = 0; i < N; i++) {
        m[0][i] = i+1;
    }
    return m;
}
Matrix MatrixCreator::diagonal_3(int N, double a1, double a2) {
    Matrix mat(N);
    mat[0][0] = a1;
    mat[0][1] = a2;
    mat[N - 1][N - 2] = a2;
    mat[N - 1][N - 1] = a1;

    for (int i = 1; i < N - 1; i++) {
        mat[i][i - 1] = a2;
        mat[i][i] = a1;
        mat[i][i + 1] = a2;
    }
    return mat;
}
Matrix MatrixCreator::diagonal_1(int N, double a1) {
    Matrix mat(N);

    for (int i = 0; i < N; i++) {
        mat[i][i] = a1;
    }
    return mat;
}