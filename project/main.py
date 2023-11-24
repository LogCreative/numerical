import numpy as np


class Solver:
    """
    求解线性方程 Ax=b
    """
    def __init__(self, A: np.array, b: np.array):
        """
        A: n x m 矩阵
        b: n x 1 矩阵
        """
        self.A = A
        self.b = b
        assert len(self.A) == len(self.b)
        self.n = len(self.A)
        self.m = len(self.A[0])

    def solve(self):
        raise NotImplementedError


class SeqGuassSolver(Solver):
    """
    顺序 Guass 消元法
    """
    def solve(self):
        A = np.copy(self.A)
        b = np.copy(self.b)

        # 消元过程
        for k in range(self.m):
            for i in range(k + 1, self.n):
                coeff = A[i][k] / A[k][k]
                for j in range(self.m):
                    A[i][j] -= A[k][j] * coeff
                b[i] -= b[k] * coeff

        # 回代过程
        x = np.empty(self.m)
        for i in reversed(range(self.n)):
            s = b[i]
            for j in range(i + 1, self.m):
                s -= A[i][j] * x[j]
            x[i] = s / A[i][i]

        return x


class PivotGuassSolver(Solver):
    """
    选主元 Guass 消元法
    """
    def solve(self):
        A = np.copy(self.A)
        b = np.copy(self.b)

        # 消元过程
        for k in range(self.m):
            # 选主元
            max_pivot = abs(A[k][k])
            max_line = k
            for p in range(k + 1, self.n):
                if abs(A[p][k]) > max_pivot:
                    max_pivot = abs(A[p][k])
                    max_line = p

            # 交换行元素
            for j in range(self.m):
                tmp = A[k][j]
                A[k][j] = A[max_line][j]
                A[max_line][j] = tmp
            tmp = b[k]
            b[k] = b[max_line]
            b[max_line] = tmp

            for i in range(k + 1, self.n):
                coeff = A[i][k] / A[k][k]
                for j in range(self.m):
                    A[i][j] -= A[k][j] * coeff
                b[i] -= b[k] * coeff

        # 回代过程
        x = np.empty(self.m)
        for i in reversed(range(self.n)):
            s = b[i]
            for j in range(i + 1, self.m):
                s -= A[i][j] * x[j]
            x[i] = s / A[i][i]

        return x


class ThomasSolver(Solver):
    """
    追赶法
    """
    def solve(self):
        A = np.copy(self.A)
        b = np.copy(self.b)

        # 检查是否是三对角矩阵以及对角占优
        is_square = self.m == self.n
        is_lower_zero = np.array_equal(np.tril(A, -2), np.zeros([self.n, self.m]))
        is_upper_zero = np.array_equal(np.triu(A, 2), np.zeros([self.n, self.m]))
        if not (is_square and is_lower_zero and is_upper_zero):
            print("不是三对角矩阵！")

        is_dominant = (abs(A[0][0]) > abs(A[0][1]) > 0)
        for i in range(self.n - 1):
            is_dominant &= (abs(A[i][i]) >= abs(A[i][i-1]) + abs(A[i][i+1]))
        is_dominant &= (abs(A[self.n-1][self.n-1]) > abs(A[self.n-1][self.n-2]) > 0)
        if not is_dominant:
            print("不是对角占优！")

        # 计算 beta
        beta = np.empty(self.n - 1)
        beta[0] = A[0][1] / A[0][0]
        for i in range(1, self.n - 1):
            beta[i] = A[i][i+1] / (A[i][i] - A[i][i-1] * beta[i-1])

        # 解 Ly=f
        y = np.empty(self.n)
        y[0] = b[0] / A[0][0]
        for i in range(1, self.n):
            y[i] = (b[i] - A[i][i-1] * y[i-1]) / (A[i][i] - A[i][i-1] * beta[i-1])

        # 解 Ux=y
        x = np.empty(self.n)
        x[self.n - 1] = y[self.n - 1]
        for i in reversed(range(self.n - 1)):
            x[i] = y[i] - beta[i] * x[i+1]

        return x


class IterativeSolver(Solver):
    """
    迭代法公共类
    """
    def __init__(self, A: np.array, b: np.array, x0 = None, maxiter = 1000, xnormthres = 1e-4):
        super().__init__(A, b)
        assert self.m == self.n
        if x0 is None:
            x0 = np.zeros(len(A))
        self.x0 = x0
        self.maxiter = maxiter
        self.xnormlist = []
        self.xnormthres = xnormthres

    def iter_manager(self, xnorm, iternum):
        if iternum > len(self.xnormlist) + 1:
            self.xnormlist.clear()
        self.xnormlist.append(xnorm)
        return xnorm >= self.xnormthres and iternum < self.maxiter


class JacobiSolver(IterativeSolver):
    """
    Jacobi 迭代法
    """
    def solve(self):
        x = np.copy(self.x0)
        A = np.copy(self.A)
        b = np.copy(self.b)
        xnorm = np.inf
        iternum = 0

        while self.iter_manager(xnorm, iternum):
            newx = np.empty(self.n)
            for i in range(self.n):
                s = 0
                for j in range(self.n):
                    if j == i:
                        continue
                    s += A[i][j] * x[j]
                newx[i] = (b[i] - s) / A[i][i]
            xnorm = np.linalg.norm(newx - x, np.inf)
            x = newx
            iternum += 1

        print("迭代次数：{}".format(iternum))
        return x


class GuassSeidelSolver(IterativeSolver):
    """
    Guass-Seidel 迭代法
    """
    def solve(self):
        x = np.copy(self.x0)
        A = np.copy(self.A)
        b = np.copy(self.b)
        xnorm = np.inf
        iternum = 0

        while self.iter_manager(xnorm, iternum):
            newx = np.empty(self.n)
            for i in range(self.n):
                s = 0
                for j in range(i):
                    s += A[i][j] * newx[j]
                for j in range(i + 1, self.n):
                    s += A[i][j] * x[j]
                newx[i] = (b[i] - s) / A[i][i]
            xnorm = np.linalg.norm(newx - x, np.inf)
            x = newx
            iternum += 1

        print("迭代次数：{}".format(iternum))
        return x


def construct_input(n):
    """
    构造输入
    """
    A = np.diag([6]*n, 0)+np.diag([8]*(n-1), -1)+np.diag([1]*(n-1), 1)
    b = np.array([7]+[15]*(n-2)+[14])
    return np.asarray(A, np.float64), np.asarray(b, np.float64)


def compare_output(x, target):
    return np.linalg.norm(x - target, np.inf)

# =================
from matplotlib import pyplot as plt


def run_solver(solver, method: str):
    for n in [10, 30, 100]:
        A, b = construct_input(n)
        target = np.array([1] * n, np.float64)
        sol = solver(A, b)
        x = sol.solve()
        norm = compare_output(x, target)
        print("{}：n={}\tinf-norm={}\n{}".format(method, n, norm, x))
        if n>50:
            solvername = solver.__name__
            diff = abs(x - target)
            plt.figure()
            plt.xlabel('$i$')
            if norm > 1000:
                plt.yscale('log')
            plt.ylabel('$|x_i-x^*_i|$')
            plt.grid(True)
            plt.plot(range(len(x)), diff)
            plt.savefig("pic/{}.pdf".format(solvername))
    return sol  # n = 100


if __name__ == '__main__':
    run_solver(SeqGuassSolver, "顺序 Guass 消元法")
    run_solver(PivotGuassSolver, "列主元 Guass 消元法")
    run_solver(ThomasSolver, "追赶法")
    jsol = run_solver(JacobiSolver, "Jacobi 迭代法")
    gssol = run_solver(GuassSeidelSolver, "Guass-Seidel 迭代法")
    plt.figure()
    plt.xlabel('$T$')
    plt.ylabel('$||x^{(k)}-x^{(k-1)}||_\infty$')
    plt.yscale('log')
    plt.grid(True)
    plt.plot(range(len(jsol.xnormlist))[1:], jsol.xnormlist[1:], label="Jacobi")
    plt.plot(range(len(gssol.xnormlist))[1:], gssol.xnormlist[1:], label="Guass-Seidel")
    plt.legend()
    plt.savefig("pic/iter.pdf")

    # cond
    for n in [10, 30, 100]:
        A, b = construct_input(n)
        print("n={}, cond(A)={}".format(n, np.linalg.cond(A)))
