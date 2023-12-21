import sympy
import math

print(sympy.simplify("5/9*(1/(2-15**(1/2)/5)+1/(2+15**(1/2)/5))+4/9"))

def f(x):
    return 1 / (x+2)

print(0.2369269*f(-0.9061793)+0.2369269*f(0.9061793)+0.4786287*f(-0.5384693)+0.4786287*f(0.5384693)+0.5688889*f(0))

a1 = 1/(-1/math.sqrt(3)+5)+1/(1/math.sqrt(3)+5)
print(a1)
a2 = 1/(-1/math.sqrt(3)+7)+1/(1/math.sqrt(3)+7)
print(a2)
a3 = 1/(-1/math.sqrt(3)+9)+1/(1/math.sqrt(3)+9)
print(a3)
a4 = 1/(-1/math.sqrt(3)+11)+1/(1/math.sqrt(3)+11)
print(a4)
print(a1+a2+a3+a4)