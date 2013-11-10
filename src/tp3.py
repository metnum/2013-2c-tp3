#coding: utf-8
# Implementación en Python del TP3 para prototipos y demas yerbas

from numpy import matrix, any
from numpy.linalg import norm, qr, solve
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

REMAIN_FACTOR = 0.95
TELEPORTATION_FACTOR = 1 - REMAIN_FACTOR


def load_data(raw_data):
    pages = int(raw_data.pop(0))
    links = int(raw_data.pop(0))

    graph = matrix(data=([[0.0] * pages] * pages))

    for link in raw_data[0: links + 1]:
        src, dest = [int(s) for s in link.split()]
        graph[src - 1, dest - 1] = 1

    #normalize by the numbe of outdegrees
    for rownum in range(graph.shape[0]):
        if norm(graph[rownum, :]) != 0:
            graph[rownum, :] = (graph[rownum, :] /
                                sum([1 if graph[rownum, i] != 0 else 0 for
                                     i in range(graph.shape[1])]))

    return graph


def has_no_out(web):
    return matrix([[1 - float(any(row)) for row in web]]).T


def v(web):
    n = web.shape[0]
    v = matrix([[1.0/ n]] * n)
    return v


def build_uniform(web):
    n, m = web.shape
    d = has_no_out(web)
    return d * v(web).T


def P_1(P, D):
    return P + D


def P_2(P, v):
    e = matrix([[1] * P.shape[0]])
    E = e.T * v.T

    c = REMAIN_FACTOR
    P_2 = c * P + (1 - c) * E
    return P_2


def pagerank_power_kamvar(P, x, epsilon):
    """
    Perform the power method for the PageRank matrix
    using Kamvar's optimization for matrix multiplication.

    P is the stochastic PageRank P' matrix
    x is the initial guess for the solution
    flow_factor 1 - TELEPORTATION_FACTOR
    """
    MAX_ITERS = 100

    k = 1
    delta = 1000000.0  # At least one iteration will be performed

    n = len(x)
    v = matrix([ 1.0/n for i in range(n) ]).T

    while delta >= epsilon and k < MAX_ITERS:
        # Optimized Ax multiplication from algorithm 1
        y = REMAIN_FACTOR * P.T * x
        w = norm(x, 1) - norm(y, 1)
        x_k = y + w * v

        delta = norm(x_k - x, 1)
        k += 1

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


def power_quad(P, x, epsilon, quad_freq):
    """
    Perform the power method for the PageRank matrix
    using Kamvar's optimization for matrix multiplication
    and quadratic optimization.

    P is the stochastic PageRank P' matrix
    x is the initial guess for the solution
    flow_factor 1 - TELEPORTATION_FACTOR
    """
    MAX_ITERS = 100

    k = 0
    delta = 1000000.0  # At least one iteration will be performed

    vec = v(web)

    x_3 = x
    x_2 = x
    x_1 = x
    x_k = x

    while delta >= epsilon and k < MAX_ITERS:
        # Optimized Ax multiplication from algorithm 1
        y = REMAIN_FACTOR * P.T * x_k
        w = norm(x_k, 1) - norm(y, 1)
        x_k = y + w * vec

        if k % quad_freq == 5:
            x_k = quad_extrapolation(x_3, x_2, x_1, x_k)
        delta = norm(x_k - x_1, 1)
        k += 1

        x_3 = x_2
        x_2 = x_1
        x_1 = x_k
        print "%s: %s" % (k, delta)
    return x_k


if __name__ == '__main__':
    filename = sys.argv[1]
    raw_data = open(filename, 'r').readlines()

    web = load_data(raw_data)

    # Create P''
    D = build_uniform(web)
    P = P_1(web, D)
    #P2 = P_2(P, v(web))

    power_quad(P, v(web), 0.00001, 15)

