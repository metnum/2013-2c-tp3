#include "assert.h"
#include <map>
#include <memory>
#include <cmath>
#include <vector>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <functional>

using namespace std;

typedef vector<double> vec;
typedef vector<vec> matrix;

typedef map<int, double> my_column;
typedef map<int, my_column> my_matrix;
typedef my_column::iterator col_iter;
typedef my_matrix::iterator matrix_iter;

template<typename T, typename ...Args>
std::unique_ptr<T> make_unique( Args&& ...args ) {
  return std::unique_ptr<T>( new T( std::forward<Args>(args)... ) );
}


double norm(vec const& v, int n=2) {
    double accum = 0.0;
    if (n == 1) {
        for(auto& val : v) {
            accum += abs(val);
        }
    } else if (n == 2) {
        for (auto& val : v) {
            accum += val * val;
        }
        accum = sqrt(accum);
    }

    return accum;
}

vec* vec_mul(double c, vec& v) {
    vec* ret = new vec(v.size());
    for(int i = 0; i < v.size(); i++) {
        (*ret)[i] = v[i] * c;
    }
    return ret;
}

vec operator -(vec const& v1, vec const& v2) {
    auto ret = vec(v1.size());
    for (int i=0; i < v1.size(); i++) {
        ret[i] = v1[i] - v2[i];
    }
    return ret;
}

void sub_inplace(vec & v1, vec const& v2) {
    for (int i=0; i < v1.size(); i++) {
        v1[i] -= v2[i];
    }
}


void sum_inplace(vec & v1, vec const& v2) {
    for (int i=0; i < v1.size(); i++) {
        v1[i] += v2[i];
    }
}

matrix operator -(matrix const& m1, matrix const& m2) {
    assert(m1.size() == m2.size());
    assert(m1[0].size() == m2[0].size());

    matrix result(m1.size(), vec(m1[0].size()));

    for (int i = 0; i < m1.size(); i++) {
        for (int j = 0; j < m1[0].size(); j++) {
            result[i][j] = m1[i][j] - m2[i][j];
        }
    }
    return result;
}


void sub_inplace(matrix& m1, matrix const& m2) {
    assert(m1.size() == m2.size());
    assert(m1[0].size() == m2[0].size());

    for (int i = 0; i < m1.size(); i++) {
        for (int j = 0; j < m1[0].size(); j++) {
            m1[i][j] -= m2[i][j];
        }
    }
}

vec operator +(vec const& v1, vec const& v2) {
    vec ret (v1.size());

    for (int i = 0; i < v1.size(); i++) {
        ret[i] = v1[i] + v2[i];
    }
    return ret;
}

vec operator /(vec& v, double c) {
    auto copied = v;
    transform(std::begin(v), std::end(v), std::begin(copied), [&](double d) { return d/c; });
    return copied;
}

void div_inplace (vec& v, double c) {
    transform(std::begin(v), std::end(v), std::begin(v), [&](double d) { return d/c; });
}

vec operator *(vec const& v1, vec const& v2) {
    auto ret = vec(v1.size());

    assert(v1.size() == v2.size());
    for (int i=0; i < v1.size(); i++) {
        ret[i] = v1[i] * v2[i];
    }
    return ret;
}

void mult_inplace(vec* v1, vec const& v2) {
    assert(v1->size() == v2.size());
    for (int i=0; i < v1->size(); i++) {
        (*v1)[i] *= v2[i];
    }
}


vec operator *(vec const& v1, double c) {
    vec ret = vec(v1);

    for (int i = 0; i < v1.size(); i++) {
        ret[i] *= c;
    }
    return ret;
}

void mult_inplace(vec& v, double c) {
    for (int i = 0; i < v.size(); i++) {
        v[i] *= c;
    }
}


vec operator *(double c, vec const& v) {
    return v * c;
}

matrix operator *(matrix const& m, double c) {
    matrix ret = matrix(m);

    for (int i = 0; i < m.size(); i++) {
        for (int j = 0; j < m[i].size(); j++) {
            ret[i][j] *= c;
        }
    }
    return ret;
}

matrix operator *(double c, matrix const& m) {
    return m * c;
}

matrix* mult_transposed(vec const& v, vec const& v_t) {
    auto res = new matrix(v.size(), vec(v_t.size()));
    for (int i=0; i < v.size(); i++) {
        for (int j=0; j < v_t.size(); j++) {
            (*res)[i][j] = v[i] * v_t[j];
        }
    }
    return res;
}

matrix* submatrix(matrix const& A) {
    matrix* result = new matrix(A.size() - 1, vec(A[0].size() - 1, 0));

    for (int i = 1; i < A.size(); i++) {
        for (int j = 1; j < A[0].size(); j++) {
            (*result)[i - 1][j - 1] = A[i][j];
        }
    }
    return result;
}

// Multiply v.T * A
vec* left_trans_multiply(vec const& v, matrix const& A) {
    vec* res = new vec(A[0].size(), 0); // Result as long as the columns in A
    for (int i=0; i < A.size(); i++) {
        for (int j=0; j < A[0].size(); j++) {
            (*res)[j] += v[i] * A[i][j];
        }
    }
    return res;
}

double left_trans_multiply(vec const& v, vec const& v2) {
    double result = 0;
    for (int i=0; i < v.size(); i++) {
        result += v[i] * v2[i];
    }
    return result;
}

matrix operator *(matrix const& a, matrix const& b) {
    int m = a[0].size();
    int n = b.size();

    auto result = matrix(m, vec(n, 0));
    for (int i = 0; i < m; i++) {
        for (int k=0; k < m; k++) {
            for (int j=0; j < n; j++) {
                result[i][j] += a[i][j] * b[k][j];
            }
        }
    }
    return result;
}

void mult_inplace(matrix& m, double d) {
    for (int i = 0; i < m.size(); i++) {
        for (int j = 0; j < m[i].size(); j++) {
            m[i][j] *= d;
        }
    }
}

void insert_block(matrix& A, matrix const& a, int row, int col) {
    for (int i = 0; i < a.size(); i++) {
        for (int j = 0; j< a[0].size(); j++) {
            A[i + row][j + col] = a[i][j];
        }
    }
}

/* Solve system of linear equations
 * assuming A is 2x2 and b is a column vector or length 2
 */
pair<double, double> solve_square_eq(matrix&  A, vec& b) {
    // System triangulated. Subtituton low
    b[1] = b[1] / A[1][1];
    double y2 = b[1];
    A[1][1] = 1;

    // Substitution hi
    b[0] -= b[1] * (A[0][1] / A[1][1]);
    double y1 = b[0] / A[0][0];

    return make_pair(y1, y2);
}


class sparse_matrix {
    private:
        my_matrix data;

    public:
        int m, n;
        double default_value;
        bool by_row;

        sparse_matrix(int m, int n, double default_value=0, bool by_row=true) {
            this->m = m;
            this->n = n;
            this->default_value = default_value;
            this->by_row = by_row;
        }

        void sum(double c) {
            this->default_value += c;
            for(auto& row: this->data) {
                for (auto& col: row.second) {
                    col.second += c;
                }
            }
        }

        void sum(sparse_matrix& m) {
            /* Default case: row-optimized left matrix, column-optimized right matrix */
            assert(!this->by_row || !m.by_row);
        }

        vec* mult(const vec& v) {
             // Multiplica esta matrix por el vector v y devuelve un puntero a
             // un nuevo vector

            vec * ret = new vec(this->n);

            assert(this->by_row); // No implementado de otra forma
            auto row_iter = this->data.begin();

            for (int i = 0; i < this->n; i++) {
                // iterador por las filas de la matriz
                double accum = 0;

                if (row_iter->first != i) {
                    if (this->default_value != 0) {
                    // La "siguiente" fila llena no es la actual
                    // Multiplico por el default
                        for(int j = 0; j < this->m; j++) {
                            accum += v[j] * this->default_value;
                        }
                    }
                } else {
                    // Esta fila esta llena
                    auto col_it = row_iter->second.begin();

                    if (this->default_value != 0) {
                        // Multiply be the empty value
                        for(int j = 0; j < this->m; j++) {
                            // Itero por las columnas de la actual fila
                            if (col_it->first == j) {
                                // La pos j esta llena
                                accum += v[j] * col_it->second;
                                col_it++;
                            } else {
                                accum += v[j] * this->default_value;
                            }
                        }
                    } else {
                        // Do not multiply by default value and a avoid a lot of ops
                        while(col_it != row_iter->second.end()) {
                            // Itero por las columnas de la actual fila
                            accum += v[col_it->first] * col_it->second;
                            col_it++;
                        }
                    }
                    row_iter++;
                }
                (*ret)[i] = accum;
            }
            return ret;
         }

        void put(int i, int j, double val) {
            assert(this->by_row);
            assert(i < n);
            assert(j < m);
            this->data[i][j] = val;
        }

        double get(int i, int j, double if_empty) {
            assert(i < n);
            assert(j < m);
            matrix_iter row = this->data.find(i);
            if (row != this->data.end()) {
                col_iter col = row->second.find(j);
                if (col != row->second.end()) {
                    return col->second;
                }
            }
            return this->default_value;
        }

        map<int, double>& get_row(int i) {
            return this->data[i];
        }

        my_matrix get_data() const{
            return this->data;
        }

};

std::ostream& operator<<(std::ostream& os, sparse_matrix const& mat) {
    for(auto& row: mat.get_data()) {
        for (auto& col: row.second) {
            os << "(" << row.first + 1 << ", " << col.first + 1 << ") -> " << col.second << endl;
        }
    }
    return os;
}

std::ostream& operator<<(std::ostream& os, vec const& v) {
    os << "[ ";
    for(auto& val: v) {
        os << val <<", ";
    }
    os << "]" << endl;
    return os;
}


