import numpy as np

def f(n):
    return n*np.sin(np.pi/n)

def Richardson():
    R = np.zeros((3,3))
    for i in range(3):
        R[i,0] = f(3*2**i)
    
    print(R[0,0].round(6), '\t')
    for i in range(1,3):
        print(R[i,0].round(6), '\t', end='')
        for j in range(1,i+1):
            R[i,j] = R[i,j-1] + (R[i,j-1]-R[i-1,j-1])/(4**j-1)
            print(R[i,j].round(6), '\t', end='')
        print()
    
    return R[2,2]

print(Richardson())
    