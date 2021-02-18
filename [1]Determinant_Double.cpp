#include <cstdio>
#include <cmath>
#include <vector>

#define  eps 1e-8

using namespace std;

//the difference of 2 lines. second line has alpha coefficient.
void difStr(vector<vector<double>> &A, int str1, int str2, int m, double alpha);

void prepCol(vector<vector<double>> &A, int col, int n, int m);//convert the column to the desired form for Gauss

int main(void) {
    int n;
    scanf("%d", &n);
    vector<vector<double>> A;
    A.resize(n, vector<double>(n ));///size n * (n + 1): another column of values.
    /*
     * Input in format:
     * a11 a12 b1
     * a21 a22 b2
     */
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            scanf("%lf", &A[i][j]);
    //to the upper triangular form.
    for (int i = 0; i < n; ++i) {
        prepCol(A, i, n, n);
    }
    ///checking determinant
    double determinant = 1;
    for (int i = 0; i < n; ++i)
        determinant *= (A[i][i]);
    printf("%d\n", (int)(round(determinant)));
}

void prepCol(vector<vector<double>> &A, int col, int n, int m) {
    double p;
    int ilast = -1;
    for (int i = col; i < n; ++i) {
        //the system is consistent - at least 1 non-zero element in the column
        if (fabs(A[i][col]) > eps && ilast == -1) {
            ilast = i;
            continue;
        }
        //we convert the column to a singular form (1,0,0....0)
        if (fabs(A[i][col]) > eps && i != ilast) {
            p = A[i][col] / A[ilast][col];
            if (i != ilast)
                difStr(A, i, ilast, m, p);
        }
    }
    //put 1 on the diagonal
    if (ilast != -1)
        A[ilast].swap(A[col]);
}

void difStr(vector<vector<double>> &A, int str1, int str2, int m, double alpha) {
    for (int i = 0; i < m; ++i) {
        A[str1][i] -= A[str2][i] * alpha;
    }
}