#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <ctime>

using namespace std;

map<set<int>, int> Minors;
vector<vector<int>> M;

/*
 * The determinant on the first n-1 lines and strings without "without".
 *For convenience, we always spread over the last row: in this case, the minors depend
 * parametrically only on columns.
 */
int DeterminantInternal(int n, set<int> strings, int without) {
    strings.erase(without);/// remove from the set of columns "without"
    if (n == 1) {///base case - 1 element
        for (auto j : strings) {
            int Minor = M[0][j];
            return Minor;
        }
    }
    if (Minors.find(strings) != Minors.end())//check to see if we have already counted this minor
        return Minors[strings];
    int ans = 0, j = 0;///j is the column number relative to the minor itself (not its number in the original matrix)
    for (auto i : strings) {

        int Minor = DeterminantInternal(n - 1, strings, i);///Minor, crossing out the i column
        int sigma = (j + n + 1) % 2;///sign of the algebraic complement
        if (sigma) Minor *= -1;
        ans += Minor * M[n - 1][i];///multiplying by algebraic complement
        ++j;
    }
    Minors[strings] = ans;///addint minor to counted
    return ans;
}

/*
 * Calculation of the determinant of a matrix using T. Laplace.
 */
int DeterminantLaplas(int n) {
    set<int> strings;///We prepare a set of lines. Initially, all lines
    for (int i = 0; i < n; ++i)
        strings.insert(i);
    return DeterminantInternal(n, strings, n);
}

void MatrixGenerator(int n) {
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            M[i][j] = rand();
}

int main() {
    int n;
    cin >> n;
    M.resize(n, vector<int>(n));
    MatrixGenerator(n);
    int t = clock();
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            cin >> M[i][j];
    cout << DeterminantLaplas(n);
    cout << '\n' << clock() - t;
    return 0;
}