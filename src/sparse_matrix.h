#include "assert.h"
#include <map>
#include <cmath>
#include <vector>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <functional>

using namespace std;

typedef vector<double> vec;
typedef vector<vec> matrix;

double norm(vec const& v, int n) {
    double accum = 0;
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

vec operator -(vec& v1, vec& v2) {
    auto ret = vec(v1.size());
    for(int i=0; i < v1.size(); i++) {
        ret[i] = v1[i] - v2[i];
    }
    return ret;
}

vec operator /(vec& v, double c) {
    auto copied = v;
    transform(std::begin(v), std::end(v), std::begin(copied), [&](double d) { return d/c; });
    return copied;
}

vec operator *(vec const& v1, vec const& v2) {
    auto ret = vec(v1.size());

    assert(v1.size() != v2.size());
    for (int i=0; i < v1.size(); i++) {
        ret[i] = v1[i] * v2[i];
    }
    return ret;
}


vec operator *(vec const& v1, double c) {
    auto ret = vec(v1.size());

    for (int i = 0; i < v1.size(); i++) {
        ret[i] = v1[i] * c;
    }
    return ret;
}

/*
matrix operator *(matrix const& m, double c) {
    auto ret = new matrix(m.size(), vec(m[0].size(), 0));

    for (int i = 0; i < m.size(); i++) {
        for (int j = 0; j < m[i].size(); j++) {
            (*ret)[i][j] += c;
        }
    }
    return ret;
} */

vec operator +(vec const& v1, vec const& v2) {
    auto ret = vec(v1.size());
    for (int i = 0; i < v1.size(); i++) {
        ret[i] = v1[i] + v2[i];
    }
    return ret;
}

class sparse_matrix {
    private:
        typedef map<int, double> my_column;
        typedef map<int, my_column> my_matrix;
        typedef my_column::iterator col_iter;
        typedef my_matrix::iterator matrix_iter;
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

        /*
        void mult_inplace(double c) {
            for (auto& row_iter : this->data) {
                for (auto& col_it : row_iter.second) {
                    col_it.second = col_it.second * c;
                }
            }
        } */

        vec* mult(const vec& v) {
             // Multiplica esta matrix por el vector v y devuelve un puntero a
             // un nuevo vector

            vec * ret = new vec(this->n);

            if (this->default_value != 0 && this->by_row) {
                // No entiendo que hace ninguno de estos dos chequeos, porque
                // chequeas que el default_value sea diferente de 0? y lo del
                // by_row ni idea :P
                auto row_iter = this->data.begin();

                for (int i = 0; i < this->n; i++) {
                    // iterador por las filas de la matriz

                    double accum = 0;

                    if (row_iter->first != i) {
                        // La "siguiente" fila llena no es la actual
                        // Multiplico por el default
                        for(int j = 0; j < this->m; j++) {
                            accum += v[j] * this->default_value;
                        }

                    } else {
                        // Esta fila esta llena
                        auto col_iter = row_iter->second.begin();

                        for(int j = 0; j < this->m; j++) {
                            // Itero por las columnas de la actual fila
                            if (col_iter->first == j) {
                                // La pos j esta llena
                                accum += v[j] * col_iter->second;
                                next(col_iter);
                            } else {
                                accum += v[j] * this->default_value;
                            }
                        }
                        next(row_iter);
                    }
                    (*ret)[i] = accum;
                }
                return ret;
            }
            assert("Multiplication method not supported");
         }

        void put(int i, int j, double val) {
            assert(this->by_row);
            this->data[i][j] = val;
        }

        double get(int i, int j, double if_empty) {
            matrix_iter row = this->data.find(i);
            if (row != this->data.end()) {
                col_iter col = row->second.find(j);
                if (col != row->second.end()) {
                    return col->second;
                }
            }
            return this->default_value;
        }

        sparse_matrix& get_column(int i);

        /*
        matrix* transponse(bool in_place=true) {
            assert(in_place);

            int i, j;
            for(auto& row: this->data) {
                for (auto& col: row.second) {
                    i = row.first;
                    j = col.first;
                    val = col.second;
                }
            }

        } */

        /*
        void row_divide(int i, double val);
        void row_sum(int i, double c);
        void row_sum(int i, sparse_matrix& row);
        void row_sub(int i, sparse_matrix& row);
        void col_sum(int j, double c);
        void col_sum(int j, sparse_matrix& col);
        void col_sub(int j, sparse_matrix& col);
        */
};

