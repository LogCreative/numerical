import numpy as np
import array_to_latex as a2l

A = np.array([[1,2,0],[2,-1,1],[0,1,3]])
B = A

Qt = np.eye(3)

for i in range(1,21):
    Q, R = np.linalg.qr(B)
    B = R @ Q
    Qt = Qt @ Q
    with np.printoptions(precision=5, suppress=True):
        print(f"{i} & \t ${a2l.to_ltx(B, frmt = '{:6.5f}', arraytype = 'pmatrix', print_out=False)}$ \\\\")

print()

with np.printoptions(precision=5, suppress=True):
    print(a2l.to_ltx(Qt, frmt = '{:6.5f}', arraytype = 'pmatrix', print_out=False))

print()
print(np.linalg.eig(A))