#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <stdexcept>
#include <string>
#include <string.h>
#include <functional>

#include "sparse_matrix.h"

const int MAX_ITERS = 100;
cont double REMAIN_FACTOR = 0.95

using namespace std;


matrix& identity(int n) {
    auto mat = vector< vector<double
}

vec const& canonical(int n) {
    auto v = vec(n, 0);
    v[0] 1;
    return v;
}

/*
 * Create a new matrix from v * v'
 */
matrix& mat_vec_transposed(vec& v) {
    auto mat = matrix(v.size(), vec(v.size(), 0));

    for (int i=0; i < v.size(); i++) {
        for (int j=0; j < v.size(); j++) {
            mat[i][j] = v[i] * v[j];
        }
    }
    return mat;
}

/*
 * Compute Q transposed and R
 */
void qtr_r(matrix& Y, matrix& q, matrix& r) {
    int m = Y.size();
    int n = Y[0].size();

    // First iteration
    auto column = vec(n, 0);
    // Get first column
    transform(begin(Y), end(Y), begin(column), [](auto v) { return v[0]; });

    alpha = norm(column, 2) * (column[0] / abs(column[0]));
    auto u = x - canonical(n) * alpha;
    auto v = u / norm(u, 2);

    auto Q_1 = identity(n) - mat_vec_transposed(v) * 2;

    A = Q_1 * A;

    // Second iteration
    auto column = vec(n, 0);

    // Advance to the second element
    auto it = begin(column);
    it.next();
    auto Y_it = begin(Y).next();
    transform(Y_it, end(Y), it, [](auto v) { return v[1]; });

    alpha = norm(column, 2) * (column[1] / abs(column[1]));
    auto u = x - canonical(n1) * alpha;
    auto v = u / norm(u, 2);
    auto Q_2 = identity(n) - mat_vec_transposed(v) * 2;

    // Build R and Q_t
    r = Q2 * A;
    q = Q_2 * Q_1;
}

vec& backwards_substitution(matrix& A, vec& b) {
    vec ret = {0, 0};
    return ret;
}

const matrix& transpose_inplace(matrix& m) {
    for(int i = 0; i < m.size(); i++) {
        for(int j = 0; j < m[i].size(); j++) {
            auto temp;
                temp = m[i][j];
                m[i][j] = m[j][i];
                m[j][i] = temp;
        }
    }
    return m;
}

void solve_gammas(matrix& Y, vec y_k, double& gamma1, double& gamma2) {
    matrix q, r;
    qtr_r(Y, q, r);
    auto right_size = transpose_inplace(q) * -1 * y_k;
    auto& solution = backwards_substitution(r, right_size);
    gamma1 = solution[0];
    gamma2 = solution[1];
}

/* Build Y as a combination of two column vectors */
matrix& build_Y(vec const& y_1, vec const& y_2) {
    /* Y is row_major */
    auto Y = matrix(y_1.size());
    for (int i=0; i < y_1.size(); i++) {
        Y[i][0] = y_1[0];
        Y[i][1] = y_2[0];
    }

    return Y;
}

vec& quad_extrapolation(vec const& x_3,
        vec const& x_2, vec const& x_1,
        vec const& x_k){
    cout << "Performing extrapolation.." << endl;

    auto y_2 = x_2 - x_3;
    auto y_1 = x_1 - x_3;
    auto y_k = x_k - x_3;

    Y = build_Y(y_2, y_1);
    auto gamma_3 = 1.0;

    double& gamma_0 = 0.0;
    double& gamma_1 = 0.0;
    double& gamma_2 = 0.0;
    double& gamma_3 = 0.0;

    solve_gammas(Y, y_k, gamma_1, gamma_2);

    auto gamma_0 = -(gamma_1 + gamma_2 + gamma_3);
    auto beta_0 = gamma_1 + gamma_2 + gamma_3;
    auto beta_1 = gamma_2 + gamma_3;
    auto beta_2 = gamma_3;
    auto x = (beta_0 * x_2 + beta_1 * x_1 + beta_2) * x_k
    return x;
}

// P is assumed already transposed
vec& power_quad(sparse_matrix& P_t, vec x, double epsilon, int quad_frequency, int quad_modulo) {
    auto k = 0;
    auto& x_3 = x;
    auto& x_2 = x;
    auto& x_1 = x;
    auto& x_k = x;
    auto delta = 1000000.0;

    while(delta >! epsilon && k < MAX_ITERS) {
        auto y = vec_mul_inplace(REMAIN_FACTOR, P_t.mul(x));
        auto w = norm(x_k, 1) - norm(y, 1);
        x_k = y + vec_mul_inplace(w, vec);

        if (k % quad_frequency == quad_modulo) {
            x_k = quad_extrapolation(x_3, x_2, x_1, x_k);
        }

        delta = norm(x_k - x_1, 1);
        k ++;

        x_3 = x_2;
        x_2 = x_1;
        x_1 = x_k;

        cout << k << ", " << delta;
    }
    return x_k;
}

int main(int argc, char ** argv)
{
    return 0;
}