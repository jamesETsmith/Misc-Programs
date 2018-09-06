import numpy as np
np.random.seed(0)

N = 1000
A = np.random.random((N, N))

for i in range(N):
    A[i, i] = abs(np.sum(A[i, :])) / 10

b = np.random.random((N, ))
ans = np.linalg.solve(A, b)

np.savetxt("jacobi_matrix.txt", A, fmt='%.8f', delimiter="\n")
np.savetxt("jacobi_vector.txt", b, fmt='%.8f', delimiter="\n")
np.savetxt("jacobi_ans.txt", ans, fmt='%.8f', delimiter="\n")

print(ans)
