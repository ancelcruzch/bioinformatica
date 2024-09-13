#include <vector>
#include <iostream>
#include <fstream>
#include <climits>
#include <cstring>

using namespace std;

// Converts base pairs and checks if they are a valid pair
int couple(char x, char y) {
    switch (x) {
        case 'A': return (y == 'U') ? -1 : 0;
        case 'U': return (y == 'A') ? -1 : 0;
        case 'G': return (y == 'C') ? -1 : 0;
        case 'C': return (y == 'G') ? -1 : 0;
    }
    return 0;
}

// Initializes the matrix with zeros
vector<vector<int>> init_matrix(int size) {
    return vector<vector<int>>(size, vector<int>(size, 0));
}

// Fills the matrix according to the RNA folding algorithm
void fill_matrix(vector<vector<int>>& nm, const string& rna) {
    int M = rna.length();
    for (int k = 1; k < M; ++k) {
        for (int i = 0; i < M - k; ++i) {
            int j = i + k;
            int left = nm[i][j - 1];
            int down = nm[i + 1][j];
            int diag = nm[i + 1][j - 1] + couple(rna[i], rna[j]);
            int rc = INT_MAX;
            for (int t = i + 1; t < j; ++t) {
                rc = min(rc, nm[i][t] + nm[t + 1][j]);
            }
            nm[i][j] = min({left, down, diag, rc});
        }
    }
}

// Traces back the minimum energy path through the matrix
void traceback(const vector<vector<int>>& nm, const string& rna, vector<pair<int, int>>& fold, int i, int j) {
    if (i < j) {
        if (i < j - 2 && nm[i][j] == nm[i + 1][j - 1] + couple(rna[i], rna[j])) {
            fold.push_back({i, j});
            traceback(nm, rna, fold, i + 1, j - 1);
        } else if (nm[i][j] == nm[i][j - 1]) {
            traceback(nm, rna, fold, i, j - 1);
        } else if (nm[i][j] == nm[i + 1][j]) {
            traceback(nm, rna, fold, i + 1, j);
        } else {
            for (int k = i + 1; k < j; ++k) {
                if (nm[i][j] == nm[i][k] + nm[k + 1][j]) {
                    traceback(nm, rna, fold, i, k);
                    traceback(nm, rna, fold, k + 1, j);
                    break;
                }
            }
        }
    }
}

// Converts the RNA and fold pairs to dot-bracket notation
string dot_write(const string& rna, const vector<pair<int, int>>& fold) {
    string dot(rna.length(), '.');
    for (auto [i, j] : fold) {
        dot[i] = '(';
        dot[j] = ')';
    }
    return dot;
}

int main() {
    string rna = "GGAAAUCC";
    auto nm = init_matrix(rna.size());
    fill_matrix(nm, rna);

    vector<pair<int, int>> fold;
    traceback(nm, rna, fold, 0, rna.size() - 1);
    string dot_bracket = dot_write(rna, fold);

    ofstream file("rna_structure_output.txt");
    file << "RNA SECUENCIA: " << rna << endl;
    file << "MATRIZ DE ENERGIA:\n";
    for (const auto& row : nm) {
        for (int num : row) {
            file << num << '\t';
        }
        file << '\n';
    }

    file << '\n';
    int n = rna.size();
    file << "Score: " << nm[0][n-1];

    file.close();

    return 0;
}
