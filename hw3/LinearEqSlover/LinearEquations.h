#ifndef LINEAREQSLOVER_LINEAREQUATIONS_H
#define LINEAREQSLOVER_LINEAREQUATIONS_H

#include "std_lib_facilities.h"
#define UNKNOWNS 101					// 最大未知数个数，x_0留空，方程组中 n 不需要连续，但是都小于100。
#define EPSILON 1e-6					// 0阈值，消除-0.0000
#define PRECISION 4

class LinearEquations {
    friend string cutDouble(double x);          // 截位函数
public:
    // 程序的输入共有N行，每一行为一个线性方程，未知数以x_n形式表示未知数。
    LinearEquations(string str);
    // “列主元高斯消去法”求解
    // pivot 参数将决定使用是否使用列主元法
    void solve(bool pivot = true, bool shorten = false);

    ~LinearEquations();
private:
    bool error = false;							// 错误的线性方程组
    int line = 1;								// 行数
    bool validUnks[UNKNOWNS] = { false };		// 有效的未知数
    double A[UNKNOWNS][UNKNOWNS] = { 0 };		// 待求矩阵	[1~maxvalidUnk][1~maxvalidUnk]
    double b[UNKNOWNS] = { 0 };					// 常数阵	[1~maxvalidUnk]
    double x[UNKNOWNS];							// 解集		[1~maxvalidUnk]
};

#endif //LINEAREQSLOVER_LINEAREQUATIONS_H
