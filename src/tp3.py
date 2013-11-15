#coding: utf-8
# Implementación en Python del TP3 para prototipos y demas yerbas

import itertools
import sys

from numpy import matrix, sign, zeros, empty
from numpy.linalg import norm, solve, qr
from scipy.sparse import lil_matrix

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
    """
    Carga una matrix con la estructura de una web dada

    El formato de archivo es:

        <n cantidad de paginas>
        <m cantidad de links>
        <l0 from> <l0 to>  link de la pagina l0.from a l0.to
        ...
        <lm from> <lm to>  link de la pagina lm.from a lm.to

    Devuelve una matriz esparsa donde la columna i indica la probabilidad de
    que un surfer en la pagina i vaya a cada pagina de la columna.

    M[i, j] = probabilidad de estando en la pagina j ir a la pagina i.

    Notar que estamos creando la matriz A de Kamvar2003 o sea, es transpuesta.
    """
    pages = int(raw_data.readline())
    links = int(raw_data.readline())

    print 'Armando la matriz esparsa...'
    #graph = zeros((pages, pages))
    graph = lil_matrix((pages, pages), dtype='d')

    row_outdegrees = zeros(pages)
    count = 0
    for i in range(links):
        # leo cada linea del archivo
        # corto en el espacio entre numeros
        # casteo a int
        # empiezo por 0 en vez de por 1
        src, dest = [int(s) - 1 for s in raw_data.readline().split()]
        graph[dest, src] = 1
        row_outdegrees[src] += 1

        # helper para ver por donde voy
        count += 1
        if count % 10000 == 0:
            print "%s of %s links loaded..." % (count, links)

    print "Normalizing out edges..."
    # get the used rows and cols
    rows, cols = graph.nonzero()
    count = 0
    for i, j in itertools.izip(rows, cols):
        graph[i, j] = 1.0 / row_outdegrees[j]
        # this can be optimized to only perform the div once per row

        # helper para ver por donde voy
        count += 1
        if count % 10000 == 0:
            print "%s of %s links normalized..." % (count, links)

    print "Matrix loaded. converting to CSR..."
    # we convert it to CSR for faster arithmetic
    graph = graph.tocsr()

    return graph, pages, links, row_outdegrees


#def D(web):
#    return matrix([[1 - float(any(row)) for row in web]]).T

def v(n):
    v = empty(n)
    v.fill(1.0 / n)
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
    # v = matrix([1.0 / n] * n).T  # we can just multiply by the const 1/n
    uniform_prob = 1.0 / n

    while delta >= epsilon and k < MAX_ITERS:
        # Optimized Ax multiplication from algorithm 1
        y = P.dot(REMAIN_FACTOR * x)
        w = norm(x, 1) - norm(y, 1)
        x_k = y + w * uniform_prob

        if criteria == 'abs':
            delta = norm(x_k - x)
        else:
            delta = norm(x_k - x) / norm(x)

        print "Iter %s, delta=%s..." % (k, delta)
        k += 1
        x = x_k

    return x_k


def solve_gammas(Y, y_k):
    # q, r = qr(Y)
    # ret = solve(r, -q.T * matrix(y_k).T).T
    # return ret.tolist()[0]
    return qr_two_iterations(Y, matrix(y_k).T).T.tolist()[0]


def quad_extrapolation(x_3, x_2, x_1, x_k):
    print "Performing interpolation..."
    y_2 = x_2 - x_3
    y_1 = x_1 - x_3
    y_k = x_k - x_3

    # Y2 = matrix([y_2[:, 0].ravel().tolist()[0], y_1[:, 0].ravel().tolist()[0]]).T
    Y = matrix((y_2, y_1)).T
    gamma_3 = 1

    (gamma_1, gamma_2) = solve_gammas(Y, y_k)

    # gamma_0 = -(gamma_1 + gamma_2 + gamma_3)
    beta_0 = gamma_1 + gamma_2 + gamma_3
    beta_1 = gamma_2 + gamma_3
    beta_2 = gamma_3
    x = beta_0 * x_2 + beta_1 * x_1 + beta_2 * x_k
    return x


def power_quad(P, x, criteria, epsilon, quad_freq=10, quad_modulo=8):
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

    # vec = v(web.shape[0])  # we can just multiply by the const 1/n
    prob_teleport = 1.0 / P.shape[0]

    x_3 = x
    x_2 = x
    x_1 = x
    x_k = x

    while delta >= epsilon and k < MAX_ITERS:
        # Optimized Ax multiplication from algorithm 1
        y = P.dot(REMAIN_FACTOR * x_k)
        w = norm(x_k, 1) - norm(y, 1)
        x_k = y + w * prob_teleport

        if k % quad_freq == quad_modulo:
            print "Performing quad extrapolation..."
            x_k = quad_extrapolation(x_3, x_2, x_1, x_k)

        if criteria == 'abs':
            delta = norm(x_k - x_1)
        else:
            delta = norm(x_k - x_1) / norm(x_1)

        x_3 = x_2
        x_2 = x_1
        x_1 = x_k
        print "Iter %s, delta=%s..." % (k, delta)
        k += 1
    return x_k


# Just two iterations of it
def qr_two_iterations(A, b):
    from ipdb import set_trace; set_trace();
    m, n = A.shape

    # First iteration
    x = A[:, 0]

    v = zeros((m, 1))
    v[0] = 1.0

    alpha = sign(x[0, 0]) * norm(x, 2)

    v = x - (alpha * v)

    v = v / norm(v, 2)

    A = A - (2 * v * (v.T * A))
    b = b - (2 * v * (v.T * b))

    # Second iteration
    x = A[1:m, 1:n][:, 0]

    v = zeros((m - 1, 1))
    v[0] = 1.0

    alpha = sign(x[0, 0]) * norm(x, 2)

    v = x - (alpha * v)

    v = v / norm(v, 2)

    A[1:m, 1:n] = A[1:m, 1:n] - (2 * v * (v.T * A[1:m, 1:n]))
    b[1:m] = b[1:m] - (2 * v * (v.T * b[1:m]))

    return solve(A[0:2, 0:2], -b[0:2])

if __name__ == '__main__':
    filename = sys.argv[1]
    raw_data = open(filename, 'r')

    print "Loading data..."
    web, pages, links, out_degrees = load_data(raw_data)

    # Create P''
    # dense = web.todense()
    # P2 = P_2(P_1(dense, out_degrees), v(len(out_degrees)))

    # print "Computing PageRank with regular Power Method..."
    # res = pagerank_power_kamvar(web, v(pages), 'rel', 0.0001)
    # print 'result:\n', res

    print
    print "Computing PageRank with regular Power Quad..."
    res = power_quad(web, v(pages), 'abs', 0.000001)
    print res / norm(res, 2)  # es necesario normalizar?

    # Uncomment to test result
    # print P2.T.dot(res) / norm(res, 2)
