#ifndef _MATRIZRALA_H
#define	_MATRIZRALA_H

#include <map>
#include <vector>
#include <iostream>

using namespace std;

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
            for(auto &kv: this->data) {
                kv.second += c;
            }
        }

        void sum(sparse_matrix& m) {
            /* Default case: row-optimized left matrix, column-optimized right matrix */
            if (!this->by_row || !m->by_row) {
                assert("Only row-optimized source and column-optimized target supported");
            }
        }

        void sub(sparse_matrix& m);
        void mult(double c);
        void mult(vector<double> vec) {
            auto ret = vector<double>(this->n);

            if (this->default_value != 0 && this->by_row) {
                auto row_iter = this->data.begin();

                for (int i = 0; i < this->n; i++) {
                    double accum = 0;

                    if (row_iter->first != i) {
                        // Empty sparse row, multiply by default value
                        for(int j = 0; j < this->m; j++) {
                            accum += vec[j] * this->default_value;
                        }

                        ret[i] = accum;
                    } else {
                        auto col_iter = this->data[i].begin();

                        for(int j = 0; j < this->m; j++) {
                            if (col_iter->first == j) {
                                accum += vec[i] * col_iter->second;
                                std::next(col_iter);
                            } else {
                                accum += vec[i] * this->default_value;
                            }
                        }

                        std::next(row_iter);
                    }
                }
                return ret;
            }

            assert("Multiplication method not supported");
        }

        sparse_matrix& mult(sparse_matrix& m2, new_by_col=true) {
            auto ret = sparse_matrix(this->m, m2->n, default_value=0, !new_by_col);

            /* Default case: row-optimized left matrix, column-optimized right matrix */
            //if (this->by_row && !m->by_row) {
                /* We have to multiply all the elements */
                if(this->default_value != 0) {

                }
            }

            return ret;
        }; // modifica m in-place

        double norm(int norm); // 1 o 2
        void put(int i, int j, double val) {
            if (!this->by_row) {
                assert("Column ordered matrix not supported");
            }

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

}

