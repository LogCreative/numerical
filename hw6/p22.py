import numpy as np

A = np.array([[5 ,5327],[5327,7277699]])
b = np.array([271.4,369321.5])
c = np.linalg.solve(A,b)
print(c)

x = [19,25,31,38,44]
y = [19.0,32.3,49.0,73.3,97.8]

def g(x):
    return c[0]+c[1]*x**2

sum = 0
for i in range(5):
    sum += (g(x[i])-y[i])**2
print(sum**(0.5))