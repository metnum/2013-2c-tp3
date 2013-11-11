#coding: utf-8
# Implementación en Python del TP3 para prototipos y demas yerbas

from numpy import matrix, any, zeros
from numpy.linalg import norm, qr, solve
from scipy.sparse import dok_matrix, lil_matrix
import sys

# Operaciones
# load_file(string filename)
# build_matrix(const double[][]& raw_data)

# sparse_matrix(int m, int n, bool by_row=true)

# Operaciones matriciales

# void set_empty_value(double c)
# void sum(double c)
# vouid sum(sparse_matrix& m)
# void sub(sparse_matrix& m)
# void mult(double c)
# void  mult(sparse_matrix& m2) // modifica m in-place
# double norm(int norm) // 1 o 2
# double put(int i, int j, double val)
# double get(int i, int j, double if_empty)
# sparse_matrix& get_column(int i);
# count_nonzero_row

# row_divide(int i, double val)
# row_sum(int i, double c)
# row_sum(int i, sparse_matrix& row)
# row_sub(int i, sparse_matrix& row)
# col_sum(int j, double c)
# col_sum(int j, sparse_matrix& col)
# col_sub(int j, sparse_matrix& col)

# ----------------------------

# Otras funciones matriz
# sparse_matrix& identity(int m, int n)

REMAIN_FACTOR = 0.85
TELEPORTATION_FACTOR = 1 - REMAIN_FACTOR


def load_data(raw_data):
    pages = int(raw_data.readline())
    links = int(raw_data.readline())

    #graph = zeros((pages, pages))
    graph = lil_matrix((pages, pages), dtype='d')

    row_outdegrees = [0] * pages
    count = 0
    for link in raw_data:
        src, dest = [int(s) for s in link.split()]
        graph[src - 1, dest - 1] = 1
        row_outdegrees[src - 1] += 1
        count += 1

        if count % 10000 == 0:
            print "%s edges loaded..." % count


    print "Normalizing out edges..."
    for rownum in range(graph.shape[0]):
        if rownum % 10000 == 0:
            print "%s pages normalized..." % rownum
        if row_outdegrees[rownum] != 0:
            for col in graph.rows[rownum]:
                graph[rownum, col] = 1.0/row_outdegrees[rownum]
        #else:
        #    graph[rownum, :] = [1.0/pages] * pages

    print "Matrix loaded. converting to CSR..."
    graph = graph.tocsr()

    return graph, row_outdegrees


#def D(web):
#    return matrix([[1 - float(any(row)) for row in web]]).T

def v(rows):
    n = rows
    v = matrix([[1.0/ n]] * n)
    return v


def D(web_outdegrees):
    n = len(web_outdegrees)
    d = matrix([[0 if row != 0 else 1 for row in web_outdegrees]]).T
    return d * v(n).T


def P_1(P, web_outdegrees):
    return P + D(web_outdegrees)


def P_2(P_1, v):
    e = matrix([[1] * P_1.shape[0]])
    E = e.T * v.T

    c = REMAIN_FACTOR
    P_2 = c * P_1 + (1 - c) * E
    return P_2


def pagerank_power_kamvar(P, x, criteria, epsilon):
    """
    Perform the power method for the PageRank matrix
    using Kamvar's optimization for matrix multiplication.

    P is the non-stochastic PageRank P matrix
    outdegrees is the vector with the count of outdegrees per
    x is the initial guess for the solution
    flow_factor 1 - TELEPORTATION_FACTOR
    """
    MAX_ITERS = 100

    k = 1
    delta = 1000000.0  # At least one iteration will be performed

    n = len(x)
    v = matrix([ 1.0/n for i in range(n) ]).T
    P_t = P.T

    while delta >= epsilon and k < MAX_ITERS:
        # Optimized Ax multiplication from algorithm 1
        y = P_t.dot(REMAIN_FACTOR * x)
        w = norm(x, 1) - norm(y, 1)
        x_k = y + w * v

        if criteria == 'abs':
            delta = norm(x_k - x)
        else:
            delta = norm(x_k - x)/ norm(x)

        print "Iter %s, delta=%s..." % (k, delta)
        k += 1
        x = x_k

    return x_k


def solve_gammas(Y, y_k):
    q, r = qr(Y)
    ret = solve(r, -q.T * y_k).T
    return ret.tolist()[0]

def quad_extrapolation(x_3, x_2, x_1, x_k):
    print "Performing interpolation..."
    y_2 = x_2 - x_3
    y_1 = x_1 - x_3
    y_k = x_k - x_3

    Y = matrix([y_2[:,0].ravel().tolist()[0], y_1[:,0].ravel().tolist()[0]]).T
    gamma_3 = 1

    (gamma_1, gamma_2) = solve_gammas(Y, y_k)

    gamma_0 = -(gamma_1 + gamma_2 + gamma_3)
    beta_0 = gamma_1 + gamma_2 + gamma_3
    beta_1 = gamma_2 + gamma_3
    beta_2 = gamma_3
    x = beta_0 * x_2 + beta_1 * x_1 + beta_2 * x_k
    return x


def power_quad(P, x, criteria, epsilon, quad_freq):
    """
    Perform the power method for the PageRank matrix
    using Kamvar's optimization for matrix multiplication
    and quadratic optimization.

    P is the stochastic PageRank P' matrix
    x is the initial guess for the solution
    flow_factor 1 - TELEPORTATION_FACTOR
    """
    MAX_ITERS = 100

    k = 1
    delta = 1000000.0  # At least one iteration will be performed

    vec = v(web.shape[0])

    x_3 = x
    x_2 = x
    x_1 = x
    x_k = x

    P_t = P.T

    while delta >= epsilon and k < MAX_ITERS:
        # Optimized Ax multiplication from algorithm 1
        y = P_t.dot(REMAIN_FACTOR * x_k)
        w = norm(x_k, 1) - norm(y, 1)
        x_k = y + w * vec

        if k % quad_freq == 5:
            print "Performing quad extrapolation..."
            x_k = quad_extrapolation(x_3, x_2, x_1, x_k)

        if criteria == 'abs':
            delta = norm(x_k - x_1)
        else:
            delta = norm(x_k - x_1)/ norm(x_1)
        k += 1

        x_3 = x_2
        x_2 = x_1
        x_1 = x_k
        print "%s: %s" % (k, delta)
    return x_k


if __name__ == '__main__':
    filename = sys.argv[1]
    raw_data = open(filename, 'r')

    print "Loading data..."
    web, out_degrees = load_data(raw_data)

    #from ipdb import set_trace; set_trace()
    # Create P''
    dense = web.todense()
    P2 = P_2(P_1(dense, out_degrees), v(len(out_degrees)))

    print "Computing PageRank with regular Power Method..."
    res = pagerank_power_kamvar(web, v(web.shape[0]), 'rel', 0.0001)
    print res / norm(res, 2)
    print "Computing PageRank with regular Power Quad..."
    res = power_quad(web, v(web.shape[0]), 'rel', 0.0001,  6)
    print res / norm(res, 2)

    # Uncomment to test result
    print P2.T.dot(res) / norm(res, 2)
