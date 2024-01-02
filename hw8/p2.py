import sympy

x = 0
y = 1

for i in range(10):
    y = 0.105 * x + 1.105 * y + 0.005
    x = x + 0.1
    print(f"${x:.1f}$\t& ${sympy.simplify(y):.4f}$\t& ${sympy.simplify(-x-1+2*sympy.exp(x)):.4f}$ \\\\")