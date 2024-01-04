import numpy as np

A = np.array([[7,3,-2],[3,4,-1],[-2,-1,3]])
x = np.array([1,1,1])

def power_method(A, x):
    vm = 0
    n = 0
    v = x
    while True:
        v = A @ v
        vmnew = np.max(v) # might be lambda
        v = v / vmnew
        n = n + 1
        print(f"{n} \t {v} \t {vmnew:.5f}")
        if (np.abs(vmnew - vm) < 1e-3):
            break
        vm = vmnew

power_method(A, x)