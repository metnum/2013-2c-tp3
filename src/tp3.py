#coding: utf-8
# Implementación en Python del TP3 para prototipos y demas yerbas

import itertools
import sys

from numpy import matrix, sign, zeros, empty
from numpy.linalg import norm, solve
from scipy.sparse import lil_matrix


REMAIN_FACTOR = 0.95
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


def pagerank_power_kamvar(P, x, criteria, epsilon, remain_factor):
    """
    Perform the power method for the PageRank matrix
    using Kamvar's optimization for matrix multiplication.

    P is the non-stochastic PageRank P matrix
    outdegrees is the vector with the count of outdegrees per
    x is the initial guess for the solution
    flow_factor 1 - TELEPORTATION_FACTOR
    """
    MAX_ITERS = 10000

    k = 1
    delta = 1000000.0  # At least one iteration will be performed

    n = len(x)
    uniform_prob = 1.0 / n

    results = [x]
    while delta >= epsilon and k < MAX_ITERS:
        # Optimized Ax multiplication from algorithm 1
        y = P.dot(remain_factor * x)
        w = norm(x, 1) - norm(y, 1)
        x_k = y + w * uniform_prob

        if criteria == 'abs':
            delta = norm(x_k - x)
        else:
            delta = norm(x_k - x) / norm(x)

        k += 1
        x = x_k
        results.append(x_k)
    for result in results[:-1]:
        print norm(abs(x_k - result), 1)

    return x_k


def solve_gammas(Y, y_k):
    return qr_two_iterations(Y, matrix(y_k).T).T.tolist()[0]


def quad_extrapolation(x_3, x_2, x_1, x_k):
    # print "Performing interpolation..."
    y_2 = x_2 - x_3
    y_1 = x_1 - x_3
    y_k = x_k - x_3

    Y = matrix((y_2, y_1)).T
    gamma_3 = 1

    (gamma_1, gamma_2) = solve_gammas(Y, y_k)

    # gamma_0 = -(gamma_1 + gamma_2 + gamma_3)
    beta_0 = gamma_1 + gamma_2 + gamma_3
    beta_1 = gamma_2 + gamma_3
    beta_2 = gamma_3
    x = beta_0 * x_2 + beta_1 * x_1 + beta_2 * x_k
    return x


def power_quad(P, x, criteria, epsilon, quad_freq, remain_factor):
    """
    Perform the power method for the PageRank matrix
    using Kamvar's optimization for matrix multiplication
    and quadratic optimization.

    P is the stochastic PageRank P' matrix
    x is the initial guess for the solution
    flow_factor 1 - TELEPORTATION_FACTOR
    """
    MAX_ITERS = 10000

    k = 1
    delta = 1000000.0  # At least one iteration will be performed

    # vec = v(web.shape[0])  # we can just multiply by the const 1/n
    prob_teleport = 1.0 / P.shape[0]

    x_3 = x
    x_2 = x
    x_1 = x
    x_k = x

    results = [x]
    while delta >= epsilon and k < MAX_ITERS:
        # Optimized Ax multiplication from algorithm 1
        y = P.dot(remain_factor * x_k)
        w = norm(x_k, 1) - norm(y, 1)
        x_k = y + w * prob_teleport

        if k in quad_freq:
            # print "Performing quad extrapolation..."
            x_k = quad_extrapolation(x_3, x_2, x_1, x_k)

        if criteria == 'abs':
            delta = norm(x_k - x_1)
        else:
            delta = norm(x_k - x_1) / norm(x_1)

        x_3 = x_2
        x_2 = x_1
        x_1 = x_k
        # print "Iter -> %s %s %s <- delta / do quad" % (k, delta, do_quad)
        k += 1
        results.append(x_k)
    for result in results[:-1]:
        print norm(abs(x_k / norm(x_k, 1) - result / norm(result, 1)), 1)
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

    print "Computing PageRank with regular Power Method... c=.9"
    res = pagerank_power_kamvar(web, v(pages), 'rel', 0.0001, 0.9)
    print 'result:\n', res

    print "Computing PageRank with regular Power Quad...c=.9, solo una vez en la 5ta iter"
    res = power_quad(web, v(pages), 'rel', 0.0001, [5], 0.9)
    print res / norm(res, 1)  # es necesario normalizar?

    print
    print "Computing PageRank with regular Power Quad...c=.9, cada 10 iter"
    res = power_quad(web, v(pages), 'rel', 0.0001, range(10, 1000, 10), 0.9)
    print res / norm(res, 1)  # es necesario normalizar?

    print
    print "Computing PageRank with regular Power Quad...c=.9, cada 5 iter"
    res = power_quad(web, v(pages), 'rel', 0.0001, range(5, 1000, 5), 0.9)
    print res / norm(res, 1)  # es necesario normalizar?

    print
    print "Computing PageRank with regular Power Quad...c=.9, cada 4 iter"
    res = power_quad(web, v(pages), 'rel', 0.0001, range(4, 1000, 4), 0.9)
    print res / norm(res, 1)  # es necesario normalizar?

    ########################### c=.95
    print
    print "Computing PageRank with regular Power Method... c=.95"
    res = pagerank_power_kamvar(web, v(pages), 'rel', 0.0001, 0.95)
    print 'result:\n', res

    print
    print "Computing PageRank with regular Power Quad...c=.95, solo una vez en la 5ta iter"
    res = power_quad(web, v(pages), 'rel', 0.0001, [5], 0.95)
    print res / norm(res, 1)  # es necesario normalizar?

    print
    print "Computing PageRank with regular Power Quad...c=.95, cada 10 iter"
    res = power_quad(web, v(pages), 'rel', 0.0001, range(10, 1000, 10), 0.95)
    print res / norm(res, 1)  # es necesario normalizar?

    print
    print "Computing PageRank with regular Power Quad...c=.95, cada 5 iter"
    res = power_quad(web, v(pages), 'rel', 0.0001, range(5, 1000, 5), 0.95)
    print res / norm(res, 1)  # es necesario normalizar?

    print
    print "Computing PageRank with regular Power Quad...c=.95, cada 4 iter"
    res = power_quad(web, v(pages), 'rel', 0.0001, range(4, 1000, 4), 0.95)
    print res / norm(res, 1)  # es necesario normalizar?

    ########################### c=.99
    print
    print "Computing PageRank with regular Power Method... c=.99"
    res = pagerank_power_kamvar(web, v(pages), 'rel', 0.0001, 0.99)
    print 'result:\n', res

    print
    print "Computing PageRank with regular Power Quad...c=.99, solo una vez en la 5ta iter"
    res = power_quad(web, v(pages), 'rel', 0.0001, [5], 0.99)
    print res / norm(res, 1)  # es necesario normalizar?

    print
    print "Computing PageRank with regular Power Quad...c=.99, cada 10 iter"
    res = power_quad(web, v(pages), 'rel', 0.0001, range(10, 1000, 10), 0.99)
    print res / norm(res, 1)  # es necesario normalizar?

    print
    print "Computing PageRank with regular Power Quad...c=.99, cada 5 iter"
    res = power_quad(web, v(pages), 'rel', 0.0001, range(5, 1000, 5), 0.99)
    print res / norm(res, 1)  # es necesario normalizar?

    print
    print "Computing PageRank with regular Power Quad...c=.99, cada 4 iter"
    res = power_quad(web, v(pages), 'rel', 0.0001, range(4, 1000, 4), 0.99)
    print res / norm(res, 1)  # es necesario normalizar?

    print
    print "Los siguientes dos casos deberian dar lo mismo y no"
    print
    print "Computing PageRank with regular Power Quad...c=.99, cada 4 iter"
    res = power_quad(web, v(pages), 'rel', 0.0001, [], 0.85)
    print res / norm(res, 1)  # es necesario normalizar?

    print
    print "Computing PageRank with regular Power Quad...c=.99, cada 4 iter"
    res = power_quad(web, v(pages), 'rel', 0.0001, [4], 0.85)
    print res / norm(res, 1)  # es necesario normalizar?


    # Uncomment to test result
    # print P2.T.dot(res) / norm(res, 2)
