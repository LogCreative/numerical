import numpy as np
import re

omega = 1.03

# SOR
DoL = np.linalg.inv(np.array([[4,0,0],[-omega,4,0],[0,-omega,4]]))
L = DoL @ np.array([[4*(1-omega),omega,0],[0,4*(1-omega),omega],[0,0,4*(1-omega)]])
f = np.array(omega * DoL @ np.array([1,4,-3]).transpose()).transpose()

x = np.array([0,0,0])
xnorm = 10000
i = 0
xfinal = np.array([0.5,1,-0.5])

while xnorm >= 5e-6:
    i = i + 1
    newx = L @ x + f
    xnorm = np.linalg.norm(xfinal-newx,np.inf)
    printx = re.sub(r'\s+',", ",str(newx)[1:-1].lstrip()).replace(".,",",")
    printxnorm = re.sub(r'e-0(\d)',r'\\times 10^{-\1}',str(np.array([xnorm]))[1:-1])
    # printx = ", ".join(map(lambda x: str(x), newx.tolist()))
    print("\\bm{{x}}^{{({})}}&=\\bm{{L}}_{{\\omega}}\\bm{{x}}^{{({})}}+\\bm{{f}}=\left({}\\right)^\\top & \\lVert\\bm{{\\epsilon}}^{{({})}}\\rVert_{{\\infty}}&={} \\\\".format(str(i), str(i-1), printx, str(i), printxnorm))
    x = newx
