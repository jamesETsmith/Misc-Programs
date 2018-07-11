#!/usr/bin/env python
'''
Light program to plot multi-column text data in a simple and quick manner.
'''

# General
import sys
import numpy as np
import matplotlib.pyplot as plt

filename = sys.argv[1]

headers = []
data = []

# Get data
with open(filename, 'r') as f:
    lines = f.readlines()

    for token in lines[0].split():
        headers.append(token)

    # Make Columns
    for col in range(len(lines[1].split())):
        data.append([])

    for i in range(1, len(lines)):
        line = lines[i].split()

        # Skip blank line
        if len(line) == 0:
            continue

        j = 0
        for token in line:
            data[j].append(float(token))
            j += 1

data = np.array(data)

for i in range(1, 3):
    data[i] -= np.mean(data[i, -10:])
    print(np.mean(data[i, -10:]))

# print(data[0])
# print(data.shape, data.size)
# Plot
plt.figure()
if data.shape[1] == data.size:
    plt.plot(np.arange(data[0].size), data[0], 'o-')
else:
    for i in range(1, data.shape[0]):
        plt.plot(data[0, :], data[i, :], label=headers[i])
    plt.xlabel(headers[0])
    plt.legend()
plt.show()
