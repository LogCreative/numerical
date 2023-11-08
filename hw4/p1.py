import numpy as np
import re

# Jacobi
B = np.array([[0,-0.4,-0.2],[0.25,0,-0.5],[-0.2,0.3,0]])
f = np.array([-2.4,5,0.3])

# Guass
# B = np.array([[ 0.   , -0.4  , -0.2  ],
#        [ 0.   , -0.1  , -0.55 ],
#        [ 0.   ,  0.05 , -0.125]])
# f = np.array([-2.4,  4.4,  2.1])

x = np.array([0,0,0])
xnorm = 10000
i = 0

while xnorm >= 1e-4:
    i = i + 1
    newx = B @ x + f
    xnorm = np.linalg.norm(newx-x,np.inf)
    printx = re.sub(r'\s+',", ",str(newx)[1:-1]).replace(".,",",")
    printxnorm = str(np.array([xnorm]))[1:-1]
    # printx = ", ".join(map(lambda x: str(x), newx.tolist()))
    print("\\bm{{x}}^{{({})}}&=\\bm{{B}}\\bm{{x}}^{{({})}}+\\bm{{f}}=\left({}\\right)^\\top & \\lVert\\bm{{\\epsilon}}^{{({})}}\\rVert_{{\\infty}}&={} \\\\".format(str(i), str(i-1), printx, str(i), printxnorm))
    x = newx