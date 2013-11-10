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

#include "sparse_matrix.h"

#define sign(X) ((X) < (0) ? (-1) : (1))

using namespace std;

// Tenemos que la matriz A
// es de m x 2 tiene una
// similar a esta:
//
//   A    x   =   b
// +--+  +-+     +-+
// |**|  |x|  =  |b|
// |**|  |y|     |c|
// |..|  +-+     +-+
// |..|
// |**|
// +--+
//
// De las cuales solo necesitamos
// hacer 0 la posicion A(2, 1)

const sparse_matrix* QR_decomposition_one_iteration(sparse_matrix& A, sparse_matrix& b)
{
    // n = size(A, 2);
    // R = A;
    // for i = 1:n
    //     x = R(:, i);

    //     e_i = zeros(n, 1);
    //     e_i(i) = 1;

    //     alpha = sign(x(i)) * norm(x, 2);

    //     u = x - (alpha * e_i);

    //     v = u / norm(u, 2);

    //     R = R - (2 * v * (v' * R));
    //     b = b - (2 * v * (v' * b)); 
    // end

    // X = R \ b;
    // return
    int m = A->m;
    int n = B->n;

    sparse_matrix& x0 = A.get_column(0);
    
    sparse_matrix& e0 = sparse_matrix.new(m, 1);
    e0.put(0, 0, 1.0);

    double alpha = sign(x0.get(0, 0)) * x0.norm(2);

    sparse_matrix& v = x0 - (x0.mult(e0));
    v = v.mult(1.0 / v.norm(2));
    
}


int main(int argc, char ** argv)
{

    

    return 0;
}
