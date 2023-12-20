// Week 1 - Problem 1:
// Romberg's Integral Method

#include<functional>
#include<iostream>
#include<iomanip>
#include<cmath>
using namespace std;


/// Romberg 积分类
class Romberg{
private:
    int         n;                          //最大细分数
    double*     h;                          //细分间距矩阵
    double      a;                          //区间起始
    double      b;                          //区间终止
    double**    R;                          //积分矩阵
    int*        BIN_P;                      //2次幂表
    int*        QUA_P;                      //4次幂表
    double      EPSILON = 1e-5;             //误差
    function<double(double)>  func;                       //函数

    // calculate parameter list h
    void calcList_h(){
        h[0] = b - a;
        for (int i = 1; i <= n; ++i)
            h[i] = h[i-1] / 2;
    }

    void calcListBIN_P(){
        BIN_P[0] = 1;
        for (int i = 1; i < n; ++i)
            BIN_P[i] = BIN_P[i-1] * 2;
    }

    void calcListQUA_P(){
        QUA_P[0] = 1;
        for (int i = 1; i <= n; ++i)
            QUA_P[i] = QUA_P[i-1] * 4;
    }

    double integral(){
        float process;
        // First, calculate R(n,0)
        R[0][0] = h[1] * (func(a) + func(b));
        for(int i = 1; i <= n; ++i){
            double sum = 0;
            for(int k = 1; k<= BIN_P[i-1]; ++k)
                sum += func(a + (2 * k - 1) * h[i]);
            R[i][0] = 0.5 * R[i-1][0] + h[i] * sum;
        }

        cout << setprecision(7) << 0 << '\t' << R[0][0] << '\n';
        // Then, calculate R(n,m)
        for(int i = 1; i<= n; ++i){
            cout << i << '\t' << R[i][0] << '\t';
            for(int j = 1; j <= i; ++j) {
                R[i][j] = R[i][j - 1] + 1.0 / (QUA_P[j] - 1) * (R[i][j - 1] - R[i - 1][j - 1]);
                cout << R[i][j] << '\t';
            }
            cout << '\n';
            if(abs(R[i][i]-R[i-1][i-1])<EPSILON)
                return R[i][i];
        }
        cout << '\n';

        return R[n][n];
    }

public:

    // constructor
    Romberg(function<double(double)> function, int maxDivsionNum){
        n = maxDivsionNum;
        h = new double[n+1];

        BIN_P = new int[n];
        calcListBIN_P();

        QUA_P = new int[n+1];
        calcListQUA_P();

        R = new double*[n+1];
        for(int i = 0; i <= n; ++i)
            R[i] = new double[n+1];

        func = function;
    }

    // distructor
    ~Romberg(){
        delete[] h;
        for(int i = 0; i <= n; ++i)
            delete[] R[i];
        delete[] R;
        delete[] BIN_P;
        delete[] QUA_P;
    }

    double Integral(int intervalStart, int intervalEnd){
        a = intervalStart;
        b = intervalEnd;
        calcList_h();
        return integral();
    }
};

double func1(double x) {
    return 2/sqrt(M_PI) * exp(-x);
}

double func2(double x) {
    return 1/x;
}

int main(){
//    function<double(double)> func1Ptr{func1};
    function<double(double)> func2Ptr{func2};

    int n;
    cout << "Please input the fineness of integral:";
    cin >> n;

    Romberg r1(func2Ptr,n);

    double a,b;
    cout << "Please input the start and end of intergral interval:";
    cin >> a >> b;

    r1.Integral(a,b);
    return 0;
}