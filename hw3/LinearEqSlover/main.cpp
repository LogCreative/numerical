#include "std_lib_facilities.h"
#include "LinearEquations.h"

int main(int argc, char* argv[]) {
    bool pivot = true, shorten = false;
    if (argc >= 2 && strcmp(argv[1], "g") == 0)
        pivot = false;
    if (argc >= 3 && strcmp(argv[2], "s") == 0)
        shorten = true;
    cout << "Current Setup: " << (pivot ? "Pivot" : "Guass") << (shorten ? ", Shortened Calculation" : "") << endl;
    cout << "One line per equation with unknown variables x_1, x_2, ...; Use Ctrl+D to finish input." << endl;

    string str;
    char ch;
    while (cin.get(ch))  // 以 EOF 作为结尾 unix: Ctrl+D
        str.push_back(ch);
    LinearEquations les(str+'\n');     // 在 linux 评测平台需要添加换行符，可能由于get(ch)的定义是不一样的

    les.solve(pivot, shorten);

    return 0;
}
