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

const sparse_matrix* QR_decomposition_two_iterations(sparse_matrix& A, sparse_matrix& b)
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
    int n = A->n;

    // First iteration
    sparse_matrix& x = A.get_column(0);
    
    sparse_matrix& e = sparse_matrix.new(m, 1);
    e.put(0, 0, 1.0);

    double alpha = sign(x.get(0, 0)) * x.norm(2);

    sparse_matrix& v = x.sub(e.mult(alpha));
    v = v.mult(1.0 / v.norm(2));

    A = A.sub(v.build_vv().mult(A).mult(2));
    b = b.sub(v.build_vv().mult(b).mult(2));
    
    // Second iteration
    x = A.get_column(1);

    e = sparse_matrix.new(m, 1);
    e.put(1, 0, 1.0);

    alpha = sign(x.get(1, 0)) * x.norm(2);

    v = x.sub(e.mult(alpha));
    v = v.mult(1.0 / v.norm(2));

    A = A.sub(v.build_vv().mult(A).mult(2));
    b = b.sub(v.build_vv().mult(b).mult(2));

    return backwards_substitution(A, b)
}

const sparse_matrix* backwards_substitution(sparse_matrix* A, sparse_matrix* b) {

}


int main(int argc, char ** argv)
{

    

    return 0;
}
