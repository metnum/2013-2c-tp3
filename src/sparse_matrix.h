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
            if (this->by_row && !m->by_row) {
            }
        }

        void sub(sparse_matrix& m);
        void mult(double c);
        void mult(sparse_matrix& m2); // modifica m in-place
        double norm(int norm); // 1 o 2
        double put(int i, int j, double val);
        double get(int i, int j, double if_empty);
        sparse_matrix& get_column(int i);

        void row_divide(int i, double val);
        void row_sum(int i, double c);
        void row_sum(int i, sparse_matrix& row);
        void row_sub(int i, sparse_matrix& row);
        void col_sum(int j, double c);
        void col_sum(int j, sparse_matrix& col);
        void col_sub(int j, sparse_matrix& col);

}

