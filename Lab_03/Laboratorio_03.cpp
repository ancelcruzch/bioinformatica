#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>

using namespace std;

string seq1 = "";
string seq2 = "";

int global_max_val = 0;

vector<vector<int> > smith_waterman(string X, string Y, int match = 1, int mismatch = -1, int gap = -2) {
    int m = X.size();
    int n = Y.size();
    vector<vector<int> > dp(m + 1, vector<int>(n + 1, 0));

    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {
            int scoreDiag = dp[i - 1][j - 1] + (X[i - 1] == Y[j - 1] ? match : mismatch);
            int scoreLeft = dp[i][j - 1] + gap;
            int scoreUp = dp[i - 1][j] + gap;
            dp[i][j] = max({scoreDiag, scoreLeft, scoreUp, 0});
        }
    }

    // Encuentra el valor máximo en la matriz dp
    int max_i = 0, max_j = 0, max_val = dp[0][0];
    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= n; j++) {
            if (dp[i][j] > max_val) {
                max_val = dp[i][j];
                max_i = i;
                max_j = j;
            }
        }
    }

    // Rastrea hacia atrás desde el valor máximo para obtener el alineamiento óptimo
    string align1 = "", align2 = "";
    while (max_i > 0 && max_j > 0) {
        if (dp[max_i][max_j] == 0) break;

        if (dp[max_i][max_j] == dp[max_i - 1][max_j - 1] + (X[max_i - 1] == Y[max_j - 1] ? match : mismatch)) {
            align1 = X[max_i - 1] + align1;
            align2 = Y[max_j - 1] + align2;
            max_i--;
            max_j--;
        } else if (dp[max_i][max_j] == dp[max_i][max_j - 1] + gap) {
            align1 = "-" + align1;
            align2 = Y[max_j - 1] + align2;
            max_j--;
        } else {
            align1 = X[max_i - 1] + align1;
            align2 = "-" + align2;
            max_i--;
        }
        cout << align1 << " " << align2 << endl;
    }

    seq1 = align1;
    seq2 = align2;
    global_max_val = max_val;

    return dp;
}

int find_index(string &X, string &Y) {
    int index;
    if(X.size()>=Y.size()){
        index = X.find(Y);
    }else{
        index = Y.find(X);
    }
    return index;
}

vector<int> find_more_index(string &cadena,string &subcadena) {
    vector<int> indices;
    int inicio = 0;
    while (true) {
        int indice = cadena.find(subcadena, inicio);
        if (indice == string::npos) {
            break;
        }
        indices.push_back(indice);
        inicio = indice + 1;
    }
    return indices;
}


int main() {

    string X = "ATGGAAGCAATATCACTGATGACTATACTACTGGTGGTAACAACAAGTAATGCAGACAAAATCTGCATCGGTCACCAATCAACAAATTCCACGGAAACTGTAGACACGCTAACAGAAACA";
    string Y = "TCAAGCGTTAGAGAAGTCATTATGTGATAAAAAAATTCAACTTGGTATCAACTTAACTAA";
    vector<vector<int> > M = smith_waterman(X, Y);

    ofstream file("output.txt");
    if(file.is_open()){
        file << "+-+-+-+-+-+-+-+-+-+SECUENCIAS+-+-+-+-+-+-+-+-+-+" << endl;
        file << "Secuencia 01: " << X << endl;
        file << "Secuencia 02: " << Y << endl;
        file << "+-+-+-+-+-+-+-+-+SMITH_WATERMAN+-+-+-+-+-+-+-+-+" << endl;
        file << "Valor Mayor: " << global_max_val << endl;
        file << "Subsecuencia en común: " << endl;
        file << seq1 << endl;
        file << seq2 << endl;
        file << "Posicion de cadenas: " << endl;
        file << find_index(X, seq1) << endl;
        file << find_index(Y, seq1) << endl;
        file << "Matriz de scores: " << endl;
        for (const auto& row : M) {
            for (int val : row) {
                file << val << " ";
            }
            file << "\n";
        }
        file << "Subsecuencias alternas: " << endl;
        vector<int> alternas = find_more_index(X, seq1);
        for (int indice : alternas) {
            file << indice << " ";
        }
        file << endl;
        file.close();
    }

    return 0;
}
