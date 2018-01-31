#!/usr/bin/env python
# Originally written by Bastien Mussard
# Modified by James Smith

import sys
import numpy as np


def get_connections(N):
    adjacency = np.zeros((N**2, N**2))
    for c in range(N**2):
        for r in range(N**2):
            # On the same row w/in one space of each other or adjacent by PBC
            if int(c / N) == int(r / N):
                # Adjacent
                if abs(r - c) == 1:
                    adjacency[r, c] = 1
                # Perdiodically adjacent
                if int(c % N) == 0 and int(r % N) == N - 1:
                    adjacency[r, c] = 1
                elif int(r % N) == 0 and int(c % N) == N - 1:
                    adjacency[r, c] = 1

            # In the same column w/in one space of each other or adjacent by PBC
            if int(c % N) == int(r % N):
                # Adjacent
                if abs(int(r / N) - int(c / N)) == 1:
                    adjacency[r, c] = 1
                # Perdiodically adjacent
                if int(c / N) == 0 and int(r / N) == N - 1:
                    adjacency[r, c] = 1
                elif int(r / N) == 0 and int(c / N) == N - 1:
                    adjacency[r, c] = 1

    # print(adjacency)
    return adjacency

def write_header(N,f):
    f.write(' &FCI NORB={:5},NELEC={:5},MS2={:5},\n'.format(N**2, N**2, N**2))
    line = '  ORBSYM='
    for i in range(N**2):
        line = line + '1,'
    f.write(line+'\n')
    f.write('  ISYM=1,\n')
    f.write(' &END\n')


def write_U_terms(N, U, f):
    for i in range(1, N**2 + 1):
        f.write('{0:5}{1:5}{2:5}{3:5}{4:13.8f}\n'.format(i, i, i, i, U))


def write_t_terms(N, adjacency, f):
    for r in range(N**2):
        for c in range(N**2):
            if adjacency[r,c] == 1:
                f.write('{0:5}{1:5}{2:5}{3:5}{4:13.8f}\n'.format(r+1, c+1, 0, 0, -1))

class Options:
    def __init__(self, filename=None, dim=None, periodic=None ):
        if filename == None: self.filename = "FCIDUMP"
        else: self.filename = filename

        if dim == None: self.dim =  "2D"
        else: self.dim = dim

        if periodic == None: self.periodic = True
        else: self.periodic = periodic


################################################################################
if len(sys.argv) == 3:
    N = int(sys.argv[1])
    adjacency = get_connections(N)
    U = float(sys.argv[2])
else:
    print('usage: ./make_fcidump_hubbard N U/|t|')
    exit()
# print N,U

options = Options()

with open(options.filename, 'w') as f:

    write_header(N,f)
    write_U_terms(N,U,f)
    write_t_terms(N,adjacency,f)

exit(0)

up_edges = []
dn_edges = []
lt_edges = []
rt_edges = []
for i in range(1, N + 1):
    up_edges.append(i)
    lt_edges.append((i - 1) * N + 1)
    rt_edges.append((i - 1) * N + N)
    dn_edges.append(N**2 - i + 1)
# print(lt_edges)
# print(up_edges)
# print(rt_edges)
# print(dn_edges)

connections = []
for i in range(1, N**2 + 1):
    if i in lt_edges:
        a = i + N - 1
    else:
        a = i - 1
    if i in up_edges:
        b = N**2 - N + 1
    else:
        b = i - N
    if i in rt_edges:
        c = i - N + 1
    else:
        c = i + 1
    if i in dn_edges:
        d = N - N**2 + i
    else:
        d = i + N
    connections.append([i, a])
    connections.append([i, b])
    connections.append([i, c])
    connections.append([i, d])
# for elt in connections:
#    print elt
print_header(N)

for i in range(1, N**2 + 1):
    print('{0:5}{1:5}{2:5}{3:5}{4:13.8f}'.format(i, i, i, i, U))
for elt in connections:
    print('{0:5}{1:5}{2:5}{3:5}{4:13.8f}'.format(elt[0], elt[1], 0, 0, -1))
