import numpy as np
import scipy

B = np.array([[0,2,1],[2,-3,1],[1,1,-5]])
P = np.array([[0,1,0],[0,0,1],[1,0,0]])

# P @ B = _ @ L @ U
_, L, U = scipy.linalg.lu(P @ B)
print(P)
print(L)
print(U)

v = np.linalg.inv(U) @ np.array([1,1,1])
vm = np.max(v)
u = v / vm

for i in range(2,20):
    y = np.linalg.solve(L, P @ u)
    v = np.linalg.solve(U, y)
    vm = np.max(v)
    u = v / vm
    print(f"{i} \t {u} \t {vm:.5f}")

print(1/vm + 6)

print(np.linalg.eig(B))