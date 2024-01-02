import sympy

a, b = sympy.symbols('a b')
x = 0
y = 0

# for i in range(10):
#     y = y + 0.1 * (a * x + b)
#     x = x + 0.1
#     print(f"${x:.1f}$\t& ${sympy.latex(sympy.simplify(y))}$\t& ${sympy.latex(sympy.simplify(0.5*a*x**2 + b*x))}$ \\\\")

for i in range(10):
    y = y + 0.1 * a * x + 0.1 * b + 0.005 * a
    x = x + 0.1
    print(f"${x:.1f}$\t& ${sympy.latex(sympy.simplify(y))}$\t& ${sympy.latex(sympy.simplify(0.5*a*x**2 + b*x))}$ \\\\")