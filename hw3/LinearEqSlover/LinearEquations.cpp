#include "LinearEquations.h"
#include "std_lib_facilities.h"

LinearEquations::LinearEquations(string str) {
    // 用户输入的合法表达式仅包含.，+，_，-符号和数字以及字母；x，
    // 每一项之间用 + 或 - 隔开，系数为整数或小数，
    // 输出小数精度设置为小数点后4位

    //对输入字符串进行分析
    int sign = 1;
    int magnit = 10;
    double dilimit = 1;
    double tmp_coeff = 0.0;
    int index = 0;
    bool indexing = false;
    char prev = '\0';
    bool right = false;

    for (auto c : str) {
        //重复特殊字符筛除
        if (((c == '.' || c == '_' || c == '+' || c == '-' || c == ' ')
             && (prev == '.' || prev == '_' || prev == '+' || prev == '-'))
            || (c == 'x' && (prev == 'x' || prev == '.' || prev == '_'))) {
            error = true; return;
        }

        if (c >= '0' && c <= '9') {
            if (indexing) index = index * 10 + c - '0';
            else {
                tmp_coeff = tmp_coeff * magnit + (c - '0') * dilimit;
                if (dilimit < 1) dilimit /= 10;
            }
        }
        else {
            switch (c) {
                case '.':
                    if (indexing) { error = true; return; }     // 不允许小数指标
                    magnit = 1; dilimit = 0.1; break;
                case 'x': tmp_coeff = sign * (tmp_coeff == 0.0 ? (prev == '0' ? tmp_coeff : 1.0) : tmp_coeff);
                    //index = 1;                                // 不允许省略 x_1
                    break;
                case '_': index = 0; indexing = true; break;
                case '+':case '-': case '=': case '\n':
                    if (tmp_coeff != 0.0) {   // 截断符
                        if (indexing) {
                            A[line][index] += tmp_coeff;
                            validUnks[index] = true;
                        }
                            //else if (index == 1)                  // 不允许省略 x_1
                            //    coeff[1] += tmp_coeff;
                        else // 常数项
                            b[line] += (right ? 1 : -1) * sign * tmp_coeff;     // 常数项分为左右
                    }
                    // 初始化参数
                    if (c == '+') sign = 1;
                    else if (c == '-') sign = -1;
                    index = 0;
                    indexing = false;
                    tmp_coeff = 0.0;
                    magnit = 10;
                    dilimit = 1;
                    if (c == '\n') {
                        if (prev == '\n') return; // 两个换行 终止
                        ++line; sign = 1; right = false; // 换行后符号、右置复位
                    }
                    if (c == '=') { right = true; sign = 1; }           // 等于后符号复位
                    break;
                default: error = true; return;
            }
        }

        prev = c;
    }
}

// 截取位数问题
// 对double小数的截位，可以通过to_string将其转换为字符串并进行截位。请自行设计相应函数。
// setprecision 为四舍五入。
// -0.0000 的问题
string cutDouble(double x) {
    //如解为0，则预取输出为0.0000
    if (abs(x-(int)x) < EPSILON) x = (int)x;
    string s = to_string(x);
    string o = "";
    int count = -1;
    for (auto c : s) {
        if (c == '.') count = 4;
        else {
            if (count == 0) return o;
            if (count > 0) --count;
        }
        o += c;
    }
}

void printDouble(double x){
    cout << setw(12) << fixed << setprecision(PRECISION) << ((abs(x-(int)x) < EPSILON) ? (int)x : x);
}

double shortenDouble(double x)
{
    double magnifier = 1.0f * pow(10, PRECISION);
    return floor(x * magnifier + 0.5f) / magnifier;
}

//考虑参杂线性相关方程组
void LinearEquations::solve(bool pivot, bool shorten) {
    // 考虑到有些未知量可能没有被使用，需要先计算有效未知数的个数
    vector<int> validIndices({ -1 });           // 有效指标，0 舍弃
    for (int i = 0; i < UNKNOWNS; ++i)
        if (validUnks[i])
            validIndices.push_back(i);
    int validUnkCnt = validIndices.size() - 1;

    // 若表达式中有不符合的输入
    // 先考虑无解的情形，行数比未知数个数少必然无解
    if (error || line - 1 < validUnkCnt) {
        cout << "error" << endl; return;
    }

    // 若有唯一解，则按下标n递增形式输出对应的解，解与解之间用一个空格隔开。
    for (int k = 1; k <= validUnkCnt; ++k) {

        // 打印
        cout << "\\begin{pmatrix}" << endl;
        for (int i = 1; i <= line - 1; ++i) {
            for (auto j = validIndices.begin() + 1; j != validIndices.end(); ++j) {
                if (shorten) A[i][*j] = shortenDouble(A[i][*j]);
                printDouble(A[i][*j]);
                cout << " & ";
            }
            cout << ' ';
            if (shorten) b[i] = shortenDouble(b[i]);
            printDouble(b[i]);
            cout << "\\\\";
            cout << endl;
        }
        cout << "\\end{pmatrix}" << endl;
        cout << endl;

        double temp;

        if (pivot) {
            // 选主元
            int maxPtr = k;
            double maxEle = abs(A[k][validIndices[k]]);
            for (int r = k + 1; r <= line - 1; ++r)                 //将所有行考虑在内
                if (abs(A[r][validIndices[k]]) > abs(maxEle)) {
                    maxPtr = r;
                    maxEle = abs(A[r][validIndices[k]]);
                }

            //交换行元素
            //for (int j = k; j <= validUnkCnt; ++j) {

            for (vector<int>::iterator j = validIndices.begin(); j != validIndices.end(); ++j) {
                temp = A[maxPtr][*j];
                A[maxPtr][*j] = A[k][*j];
                A[k][*j] = temp;
            }
            temp = b[maxPtr];
            b[maxPtr] = b[k];
            b[k] = temp;

            if (!A[k][validIndices[k]]) {
                cout << "error" << endl;
                return;
            }
        }

        for (int i = k + 1; i <= line - 1; ++i) {               //将所有行考虑在内
            temp = A[i][validIndices[k]] / A[k][validIndices[k]];
            for (vector<int>::iterator j = validIndices.begin() + k; j != validIndices.end(); ++j)
                A[i][*j] = A[i][*j] - temp * A[k][*j];
            b[i] = b[i] - temp * b[k];
        }
    }

    // 回代过程
    // end 得到数组的最后一个单元+1的指针
    x[validIndices[validUnkCnt]] = b[validUnkCnt] / A[validUnkCnt][validIndices[validUnkCnt]];
    for (int k = validUnkCnt - 1; k >= 1; --k) {
        double S = b[k];
        //for (int j = k + 1; j <= validUnkCnt; ++j)
        for (vector<int>::iterator j = validIndices.begin() + k + 1; j != validIndices.end(); ++j)
            S -= A[k][*j] * x[*j];
        x[validIndices[k]] = S / A[k][validIndices[k]];
        // Linux 平台，Windows 平台，同时考虑inf和nan
        if (isinf(x[validIndices[k]])||isnan(x[validIndices[k]])) {
            cout << "error" << endl; return;        // 无解
        }
    }

    cout << "\\begin{align*}" << endl;
    for (int k = 1; k <= validUnkCnt; ++k) {
        cout << "x_" << k << "&=";
        printDouble(x[validIndices[k]]);
        if (k < validUnkCnt) cout << " & ";
    }
    cout << endl << "\\end{align*}" << endl;

}

LinearEquations::~LinearEquations() = default;