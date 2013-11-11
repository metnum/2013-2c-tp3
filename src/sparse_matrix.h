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

double norm(vec v, int n) {
    double accum = 0;
    if (n == 1) {
        for(double& val : v) {
            accum += abs(val);
        }
    } else if (n == 2) {
        for (double& val : v) {
            accum += val * val;
        }
        accum = sqrt(accum);
    }

    return accum;
}

vec& vec_mul_inplace(double c, vec& v) {
    for(auto& val: v) {
        val *= c;
    }
    return v;
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
    transform(begin(v), end(v), begin(copied), [&](double d) { return d/c; });
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

class sparse_matrix {
    private:
        map<int, map <int, double>> data;

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

        void sub_inplace(sparse_matrix& m) {
        }

        void mult_inplace(double c) {
            for(auto& row_iter: begin(this->data)) {
                for(auto& col_iter: begin(row_iter.begin) {
                    col_iter->second = col_iter->second * c;
                }
            }
        }

        sparse_matrix& mult(double c) {
            auto& ret = sparse_matrix(this->m, this->n, this->default_value *= c, this->by_row);
            for(auto& row_iter: this->data.begin()) {
                for(auto col_iter: row_iter.begin()) {
                    ret->put(row_iter->first, col_iter->first,
                            col_iter->second * c);
                }
            }
            return sparse_matrix;
        }

        vec& mult(const vec& v) {
            auto ret = vec(this->n);

            if (this->default_value != 0 && this->by_row) {
                auto row_iter = this->data.begin();

                for (int i = 0; i < this->n; i++) {
                    double accum = 0;

                    if (row_iter->first != i) {
                        // Empty sparse row, multiply by default value
                        for(int j = 0; j < this->m; j++) {
                            accum += v[j] * this->default_value;
                        }

                        ret[i] = accum;
                    } else {
                        auto col_iter = this->data[i].begin();

                        for(int j = 0; j < this->m; j++) {
                            if (col_iter->first == j) {
                                accum += v[i] * col_iter->second;
                                std::next(col_iter);
                            } else {
                                accum += v[i] * this->default_value;
                            }
                        }
                        next(row_iter);
                    }
                }
                return ret;
            }

            assert("Multiplication method not supported");
        }

        sparse_matrix& mult(sparse_matrix& m2, bool new_by_col=true) {
            auto ret = sparse_matrix(this->m, m2->n, default_value=0, !new_by_col);

            /* Default case: row-optimized left matrix, column-optimized right matrix */
            //if (this->by_row && !m->by_row) {
                /* We have to multiply all the elements */
            assert(0);
            if(this->default_value != 0) {
            }
            return ret;
        }; // modifica m in-place

        double norm(int norm) {
        }

        void put(int i, int j, double val) {
            assert(!this->by_row);
            this->data[i][j] = val;
        }

        double get(int i, int j, double if_empty) {
            if (auto row = this->data.at(i)) {
                if (auto col = row->at(j)) {
                    return col;
                }
            }
            return this->default_value;
        }

        sparse_matrix& get_column(int i);

        void row_divide(int i, double val);
        void row_sum(int i, double c);
        void row_sum(int i, sparse_matrix& row);
        void row_sub(int i, sparse_matrix& row);
        void col_sum(int j, double c);
        void col_sum(int j, sparse_matrix& col);
        void col_sub(int j, sparse_matrix& col);
};

