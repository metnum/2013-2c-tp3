#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <memory>
#include <vector>
#include <string>
#include <functional>

#include "sparse_matrix.h"

using namespace std;

const int MAX_ITERS = 100;
const double REMAIN_FACTOR = 0.85;

/*
 * Compute Q transposed and R
 * Solve least squares problem
 */
std::pair<double, double> qr_two_iterations(matrix& Y, vec& b) {
    int m = Y.size();
    int n = Y[0].size();

    // First iteration
    auto x = make_unique<vec> (vec(m));
    // Get first column
    for (int i=0; i < m; i++) {
        (*x)[i] = Y[i][0];
    }

    // V is the first canonical vector
    unique_ptr<vec> v (new vec(m, 0));
    (*v)[0] = 1;

    auto alpha = ((*x)[0] / abs((*x)[0])) * norm((*x), 2);
    // v = x - (v * alpha);
    mult_inplace(*v, -alpha);
    sum_inplace(*v, *x);
    // v = v / ||v||
    div_inplace(*v, norm(*v, 2));

    matrix & A = Y;
    vec& r = b;

    //cout << "v.T * A" << endl;
    unique_ptr<vec> vt_A (left_trans_multiply(*v, A));
    //cout << "v.T * b" << endl;
    double vt_B = left_trans_multiply(*v, r);

    //cout << "v * (v.T * A)" << endl;
    unique_ptr<matrix> vvtA (mult_transposed(*v, *vt_A));
    //cout << "2 * v * (v.T * b)" << endl;
    unique_ptr<vec> vvtB = make_unique<vec> ((*v) * (2 * vt_B));

    //cout << "2 * (v.T * A)" << endl;
    mult_inplace(*vvtA, 2);
    //cout << "A - (2 * (v.T * A))" << endl;
    sub_inplace(A, *vvtA);

    // mult_inplace(*vvtB, 2);
    unique_ptr<vec> vvtBd = move(vvtB);
    sub_inplace(r, *vvtBd);
    //unique_ptr<vec> r_res = make_unique<vec> (r - *vvtBd);
    vec& r_res = r;
    // Second iteration
    //cout << "Allocating x2" << endl;
    unique_ptr<vec> x_2 (new vec(m - 1, 0));

    //cout << "x2 = A[1:m, 1:n]" << endl;
    for (int i=0; i < A.size() - 1; i++) {
        (*x_2)[i] = A[i + 1][1];
    }

    //cout << "Allocating v2" << endl;
    unique_ptr<vec> v_2 (new vec(m - 1, 0));
    (*v_2)[0] = 1;

    auto alpha_2 = ((*x_2)[0] / abs((*x_2)[0])) * norm((*x_2), 2);

    //cout << "v2 = v2 * alpha" << endl;
    mult_inplace(*v_2, -alpha_2);

    // v2 = x - (v2 * alpha);
    //cout << "v2 = v2 - (v2 * alpha)" << endl;
    sum_inplace(*v_2, *x_2);

    // v = v / ||v||
    //cout << "v2 = v2 / ||v2||" << endl;
    div_inplace(*v_2, norm(*v_2, 2));

    //cout << "Getting submatrix of A" << endl;
    unique_ptr<matrix> A_2 (submatrix(A));

    // Generate orthogonal stuff
    //cout << "Generating sub-vector or b" << endl;
    vec r_res_sub = vec(r_res.begin() + 1, r_res.end());

    //cout << "v2.T * A2" << endl;
    unique_ptr<vec> vt_A2 (left_trans_multiply(*v_2, *A_2));

    //cout << "v2 * b[1:m]" << endl;
    double vt_B2 = left_trans_multiply(*v_2, r_res_sub);

    //cout << "v2 * (v2.T * A2)" << endl;;
    unique_ptr<matrix> vvtA2 (mult_transposed(*v_2, *vt_A2));
    //cout << "v2 * (v2.T * B2)" << endl;;
    vec vvtB2 = vt_B2 * (*v_2);

    //cout << "2 * v2 * (v2.T * A2)" << endl;
    mult_inplace(*vvtA2, 2);

    //cout << "A2 = A2 - (2 * v2 * (v2.T * A2)" << endl;
    sub_inplace(*A_2, *vvtA2);

    //cout << "b2 = b2 - (2 * v2 * (v2.T * b2)"  << endl;
    mult_inplace(vvtB2, 2);
    sub_inplace(r_res_sub, vvtB2);
    vec& r_2 = r_res_sub;
    // Insert slice into original matrix

    //cout << "Inserting b2 into b" << endl;
    for(int i=1; i < r.size(); i++) {
        r_res[i] = r_2[i-1];
    }

    //cout << "Inserting A2 into A" << endl;
    insert_block(A, *A_2, 1, 1);

    //cout << "b = b * -1" << endl;
    mult_inplace(r_res, -1);
    vec& inverted_r = r_res;

    //cout << "Least squares" << endl;
    return solve_square_eq(A, inverted_r);
}

/* Build Y as a combination of two column vectors */
matrix* build_Y(vec const& y_1, vec const& y_2) {
    // Y is row_major
    auto Y = new matrix(y_1.size(), vec(2, 0));
    for (int i=0; i < y_1.size(); i++) {
        (*Y)[i][0] = y_1[i];
        (*Y)[i][1] = y_2[i];
    }

    return Y;
}

unique_ptr<vec> quad_extrapolation(vec const& x_3,
        vec const& x_2, vec const& x_1, vec const& x_k){
    vec y_2 = x_2 - x_3;
    vec y_1 = x_1 - x_3;
    vec y_k = x_k - x_3;

    unique_ptr<matrix> Y (build_Y(y_2, y_1));

    auto gamma_3 = 1.0;
    double gamma_1;
    double gamma_2;

    auto gammas = qr_two_iterations(*Y, y_k);
    gamma_1 = gammas.first;
    gamma_2 = gammas.second;

    auto gamma_0 = -(gamma_1 + gamma_2 + gamma_3);
    auto beta_0 = gamma_1 + gamma_2 + gamma_3;
    auto beta_1 = gamma_2 + gamma_3;
    auto beta_2 = gamma_3;

    auto b0 = beta_0 * x_2;
    auto b1 = beta_1 * x_1;
    auto b2 = beta_2 * x_k;
    unique_ptr<vec> x (new vec(b0 + b1 + b2));
    return x;
}

vec* build_uniform(int size) {
    auto value = 1.0 / size;
    return new vec(size, value);
}

// P is assumed already transposed
vec power_quad(sparse_matrix& P_t, double epsilon, int quad_frequency=10, int quad_modulo=8) { //, int quad_frequency, int quad_modulo) {
    auto k = 1;
    shared_ptr<vec> x (build_uniform(P_t.m));
    shared_ptr<vec> x_3 = x;
    shared_ptr<vec> x_2 = x;
    shared_ptr<vec> x_1 = x;
    shared_ptr<vec> x_k = x;
    vec* v = build_uniform(x->size());

    double delta = 1000000.0;

    while(delta >= epsilon && k < MAX_ITERS) {
        unique_ptr<vec> mult (P_t.mult(*x_k));
        unique_ptr<vec> y (vec_mul(REMAIN_FACTOR, *mult));

        double w = norm(*x_k, 1) - norm(*y, 1);
        unique_ptr<vec> w_v (vec_mul(w, *v));
        sum_inplace(*y, *w_v);
        x_k = move(y);
        //cout << "Res: " << *x_k;
        //cout << "X previo: " << *x_1;

        if (k % quad_frequency == quad_modulo) {
            cout << "Performing extrapolation..." << endl;
            x_k = quad_extrapolation(*x_3, *x_2, *x_1, *x_k);
        }

        // TODO: agregar criterior relativo de detencion como el python
        delta = norm(*x_k - *x_1, 2);

        x_3 = x_2;
        x_2 = x_1;
        x_1 = x_k;

        cout << "Iter " << k << ", delta=" << delta << endl;
        k++;
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
    unique_ptr<sparse_matrix> matrix  (load_matrix(filename));

    // cout << *matrix << endl;
    // vec* initial = build_uniform(matrix->n);
    auto solution = power_quad(*matrix, 0.000001);

    //cout << "Solution: " << solution << endl;
     solution = solution / norm(solution);
     cout << solution[0] << ", " << solution[1] << ", " << solution[2] << "...." << endl;
     cout << solution[solution.size() - 4] << ", " << solution[solution.size() - 3] << ", " << solution[solution.size() - 2] << ", " << solution[solution.size() - 1] << endl;
    return 0;
}
