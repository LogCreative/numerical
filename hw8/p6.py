
h = 0.2

# def f(x,y):
#     return x + y

def f(x,y):
    return 3 * y / (1 + x)

def RK4(x, y):
    K1 = f(x, y)
    K2 = f(x + 0.5 * h, y + 0.5 * h * K1)
    K3 = f(x + 0.5 * h, y + 0.5 * h * K2)
    K4 = f(x + h, y + h * K3)
    return y + h * (K1 + 2 * K2 + 2 * K3 + K4) / 6

x = 0
y = 1
for i in range(5):
    y = RK4(x, y)
    x = x + h
    print(f"{x:.1f} & {y:.4f} \\\\")