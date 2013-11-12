#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <memory>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <stdexcept>
#include <string>
#include <string.h>
#include <functional>

#include "sparse_matrix.h"

using namespace std;

const int MAX_ITERS = 100;
const double REMAIN_FACTOR = 0.95;


/*
matrix* identity(int n) {
    matrix* m = new matrix(n, vec(n, 0));
    for (int i = 0; i < n; i++) {
        (*m)[i][i] = 1;
    }

    return m;
}
*/


// Canonical vector with first position set to 1.
/*
vec* canonical(int n) {
    vec* v = new vec(n, 0);
    (*v)[0] = 1;
    return v;
}
*/

/*
 * Create a new matrix from v * v'
 */

/*
matrix* mat_vec_transposed(vec& v) {
    auto mat = new matrix(v.size(), vec(v.size(), 0));

    for (int i=0; i < v.size(); i++) {
        for (int j=0; j < v.size(); j++) {
            (*mat)[i][j] = v[i] * v[j];
        }
    }
    return mat;
}
*/

/*
 * Compute Q transposed and R
 */
/*
void qtr_r(matrix& Y, matrix& q, matrix& r) {
    int m = Y.size();
    int n = Y[0].size();

    // First iteration
    auto x = vec(n, 0);
    // Get first column
    transform(begin(Y), end(Y), begin(x), [](vector<double> v) { return v[0]; });

    auto alpha = norm(x, 2) * (x[0] / abs(x[0]));
    auto alpha_canonical = (*canonical(n)) * alpha;
    auto u = x - alpha_canonical;
    auto v = u / norm(u, 2);

    auto Q_1 = (*identity(n)) - (*mat_vec_transposed(v)*) * 2.0;
    auto A = Q_1 * A;

    // Second iteration
    auto x = vec(n, 0);

    // Advance to the second element
    auto it = begin(x );
    it.next();
    auto Y_it = begin(Y).next();
    transform(Y_it, end(Y), it, [](auto v) { return v[1]; });

    alpha = norm(x , 2) * (x [1] / abs(x [1]));
    auto u = x - canonical(n1) * alpha;
    auto v = u / norm(u, 2);
    auto Q_2 = *identity(n) - mat_vec_transposed(v) * 2;

    // Build R and Q_t
    r = Q2 * A;
    q = Q_2 * Q_1;
}
*/

/*
vec& backwards_substitution(matrix& A, vec& b) {
    vec ret = {0, 0};
    return ret;
}
*/

/*
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
*/

/*
void solve_gammas(matrix& Y, vec y_k, double& gamma1, double& gamma2) {
    matrix q, r;
    qtr_r(Y, q, r);
    auto right_size = transpose_inplace(q) * -1 * y_k;
    auto& solution = backwards_substitution(r, right_size);
    gamma1 = solution[0];
    gamma2 = solution[1];
}
*/

/* Build Y as a combination of two column vectors */
/*
matrix& build_Y(vec const& y_1, vec const& y_2) {
    // Y is row_major
    auto Y = matrix(y_1.size());
    for (int i=0; i < y_1.size(); i++) {
        Y[i][0] = y_1[0];
        Y[i][1] = y_2[0];
    }

    return Y;
}
*/

/*
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
*/

vec* build_uniform(int size) {
    auto value = 1.0 / size;
    return new vec(size, value);
}

// P is assumed already transposed
vec power_quad(sparse_matrix& P_t, double epsilon) { //, int quad_frequency, int quad_modulo) {
    auto k = 0;
    shared_ptr<vec> x (build_uniform(P_t.m));
    shared_ptr<vec> x_3 = x;
    shared_ptr<vec> x_2 = x;
    shared_ptr<vec> x_1 = x;
    shared_ptr<vec> x_k = x;
    vec* v = build_uniform(x->size());

    double delta = 1000000.0;

    while(delta >= epsilon && k < MAX_ITERS) {
        cout << "Multiplicando P_t por x_k por const" << endl;

        unique_ptr<vec> mult (P_t.mult(*x_k));
        unique_ptr<vec> y (vec_mul(REMAIN_FACTOR, *mult));

        double w = norm(*x_k, 1) - norm(*y, 1);
        unique_ptr<vec> w_v (vec_mul(w, *v));
        unique_ptr<vec> res (new vec(*y + *w_v));
        x_k = std::move(res);
        //cout << "Res: " << *x_k;
        //cout << "X previo: " << *x_1;

        /*
        if (k % quad_frequency == quad_modulo) {
            x_k = quad_extrapolation(x_3, x_2, x_1, x_k);
        }
        */

        // TODO: agregar criterior relativo de detencion como el python
        delta = norm(*x_k - *x_1, 2);
        k++;

        x_3 = x_2;
        x_2 = x_1;
        x_1 = x_k;

        cout << "Iter " << k << ", delta=" << delta << endl;
    }

    delete v;

    return (*x_k) / norm(*x_k);
}

sparse_matrix* load_matrix(string filename) {
    ifstream file(filename);

    if (!file.is_open()) {
        cout << "No se pudo abrir el archivo " << filename << endl;
        exit(1);
    }

    double n;
    double links;

    file >> n;
    file >> links;

    cout << "Loading data..." << endl;
    auto matrix = new sparse_matrix(n, n);

    int i, j;

    auto colnum_count = vector<int>(n, 0);

    for (auto k = 0; k < links; k++) {
        file >> i;
        file >> j;

        // Matrix is built with transposed indices
        matrix->put(j - 1, i - 1, 1);
        colnum_count[i -1]++;
    }

    for (auto row: matrix->get_data()) {
        for (auto col: row.second) {
            matrix->put(row.first, col.first, 1.0 / colnum_count[col.first]);
        }
    }

    return matrix;
}

int main(int argc, char ** argv) {
    if (argc < 2) {
        cout << "Usage: tp3 [matrix_filename]" << endl;
        exit(1);
    }

    string filename = argv[1];
    auto matrix = load_matrix(filename);

    cout << *matrix << endl;
    // vec* initial = build_uniform(matrix->n);
    auto solution = power_quad(*matrix, 0.0001);

    /* Test sparse mult */
    /*
    auto& matrix = sparse_matrix(
    */

    cout << "Solution: " << solution << endl;
    // cout << solution[0] << ", " << solution[1] << ", " << solution[2] << "...." << solution[solution.size() - 2] << ", " << solution[solution.size() - 1] << endl;
    //free matrix;
    delete matrix;
    return 0;
}
