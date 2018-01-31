#!/usr/bin/env python
# Originally written by Bastien Mussard
# Modified by James Smith

import sys
import numpy as np

################################################################################

def get_connections(options, N):
    adjacency = np.zeros((N**2, N**2))
    for c in range(N**2):
        for r in range(N**2):
            # On the same row w/in one space of each other or adjacent by PBC
            if int(c / N) == int(r / N):
                # Adjacent
                if abs(r - c) == 1:
                    adjacency[r, c] = 1
                if options.periodic:
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
                if options.periodic:
                    # Perdiodically adjacent
                    if int(c / N) == 0 and int(r / N) == N - 1:
                        adjacency[r, c] = 1
                    elif int(r / N) == 0 and int(c / N) == N - 1:
                        adjacency[r, c] = 1

    # print(adjacency)
    return adjacency

def write_header(options, N,f):
    f.write(' &FCI NORB={:5},NELEC={:5},MS2={:5},\n'.format(N**2, N**2, N**2))
    line = '  ORBSYM='
    for i in range(N**2):
        line = line + '1,'
    f.write(line+'\n')
    f.write('  ISYM=1,\n')
    f.write(' &END\n')


def write_U_terms(options, N, U, f):
    for i in range(1, N**2 + 1):
        f.write('{0:5}{1:5}{2:5}{3:5}{4:13.8f}\n'.format(i, i, i, i, U))


def write_t_terms(options, N, adjacency, f):
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
if __name__ == "__main__":

    if len(sys.argv) == 3:
        N = int(sys.argv[1])
        U = float(sys.argv[2])
    else:
        print('usage: ./make_fcidump_hubbard N U/|t|')
        exit()

    options = Options()
    adjacency = get_connections(options, N)
    with open(options.filename, 'w') as f:

        write_header(options,N,f,)
        write_U_terms(options,N,U,f,)
        write_t_terms(options,N,adjacency,f,)
