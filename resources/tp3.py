#coding: utf-8
# Implementación en Python del TP3 para prototipos y demas yerbas

from numpy import matrix, any
import sys


TELEPORTATION_FACTOR = 0.95


def load_data(raw_data):
    pages = int(raw_data.pop(0))
    links = int(raw_data.pop(0))

    graph = matrix(data=([[0.0] * pages] * pages))

    for link in raw_data[0: links + 1]:
        src, dest = [int(s) for s in link.split(' ')]
        graph[src - 1, dest - 1] = 1

    return graph


def has_no_out(web):
    return [1 - float(any(row)) for row in web]


def v(web):
    n = web.shape[0]
    return matrix([[1.0/ n]] * n)


def build_uniform(web):
    n, m = web.shape
    d = matrix([has_no_out(web)]).T
    return d * v(web).T


def P_1(P, D):
    return P + D


def P_2(P, v):
    e = matrix([[1] * P.shape[0]])
    E = e.T * v.T

    c = TELEPORTATION_FACTOR
    P_2 = c * P + (1 - c) * E
    return P_2


if __name__ == '__main__':
    filename = sys.argv[1]
    raw_data = open(filename, 'r').readlines()

    web = load_data(raw_data)

    # Create P''
    D = build_uniform(web)
    P = P_1(web, D)
    P2 = P_2(P, v(web))
    print P2

    # Redefine P
    P = P2.T
    # Solve
    initial = matrix([[1.0 / P.shape[0]] * P.shape[0]])

    # Assume initial
    y = TELEPORTATION_FACTOR * P * initial.T


