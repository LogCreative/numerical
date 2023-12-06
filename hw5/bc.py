import sympy

x = [0,1,2,3]
y = [0,0,0,0]
m = ['1','-4/15','1/15','0'] #(1)
# m = ['-13/45','7/90','-1/45','1/90'] #(2)

def beta(i):
    return f'(x-{x[i]})*((x-{x[i+1]})**2)'

def beta1(i):
    return f'(x-{x[i+1]})*((x-{x[i]})**2)'

def S(i):
    return sympy.latex(sympy.expand(f'{m[i]}*{beta(i)}+{m[i+1]}*{beta1(i)}'))

print(S(0))
print(S(1))
print(S(2))