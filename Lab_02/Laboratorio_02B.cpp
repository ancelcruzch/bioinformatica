/*
 ESTUDIANTE: ANCEL ALAIN FERNANDO, CRUZ CHAIÑA
 CURSO: BIOINFORMATICA GRUPO-A
 TEMA: LABORATORIO_02
 */

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <cctype>
#include <chrono>

using namespace std;

// Pregunta 03
// Implementar el alineamiento global.
void printMatrix(vector<vector<int> > &matrix) {
    for (auto& row : matrix) {
        for (auto& elem : row) {
            cout << elem << " ";
        }
        cout << "\n";
    }
}

void convertirAMayusculas(string& str) {
    for (char& c : str) {
        c = toupper(c);
    }
}

// Función para hacer el backtracking y encontrar todos los alineamientos óptimos.
void traceback(vector<vector<int>>& matrix, string& seq1, string& seq2, int i, int j, string& aligned1, string& aligned2, int match, int mismatch, int gap, vector<pair<string, string>>& alignments) {
    if (i == 0 && j == 0) {
        reverse(aligned1.begin(), aligned1.end());
        reverse(aligned2.begin(), aligned2.end());
        alignments.push_back({aligned1, aligned2});
        reverse(aligned1.begin(), aligned1.end());
        reverse(aligned2.begin(), aligned2.end());
        return ;
    }

    if (i > 0 && matrix[i][j] == matrix[i-1][j] + gap) {
        aligned1.push_back(seq1[i-1]);
        aligned2.push_back('-');
        traceback(matrix, seq1, seq2, i-1, j, aligned1, aligned2, match, mismatch, gap, alignments);
        aligned1.pop_back();
        aligned2.pop_back();
    }
    if (j > 0 && matrix[i][j] == matrix[i][j-1] + gap) {
        aligned1.push_back('-');
        aligned2.push_back(seq2[j-1]);
        traceback(matrix, seq1, seq2, i, j-1, aligned1, aligned2, match, mismatch, gap, alignments);
        aligned1.pop_back();
        aligned2.pop_back();
    }
    if (i > 0 && j > 0 && matrix[i][j] == matrix[i-1][j-1] + (seq1[i-1] == seq2[j-1] ? match : mismatch)) {
        aligned1.push_back(seq1[i-1]);
        aligned2.push_back(seq2[j-1]);
        traceback(matrix, seq1, seq2, i-1, j-1, aligned1, aligned2, match, mismatch, gap, alignments);
        aligned1.pop_back();
        aligned2.pop_back();
    }
}

void global_seq(string seq1, string seq2) {
    int match = 1;
    int mismatch = -1;
    int gap = -2;
    int n = seq1.size();
    int m = seq2.size();
    vector<vector<int>> matrix(n + 1, vector<int>(m + 1, 0));
    for (int i = 0; i <= n; ++i) {
        matrix[i][0] = i * gap;
    }
    for (int j = 0; j <= m; ++j) {
        matrix[0][j] = j * gap;
    }

    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= m; ++j) {
            int diagonal = matrix[i - 1][j - 1] + (seq1[i - 1] == seq2[j - 1] ? match : mismatch);
            int up = matrix[i - 1][j] + gap;
            int left = matrix[i][j - 1] + gap;
            matrix[i][j] = max({diagonal, up, left});
        }
    }
    cout << "Matriz de alineamiento:\n";
    printMatrix(matrix);

    vector<pair<string, string>> alignments;
    string aligned1, aligned2;
    traceback(matrix, seq1, seq2, n, m, aligned1, aligned2, match, mismatch, gap, alignments);
    cout << "Alineamientos óptimos:\n";
    for (const auto& alignment : alignments) {
        cout << alignment.first << "\n";
        cout << alignment.second << "\n";
        cout << "\n";
    }

}

int main(){
    auto start = std::chrono::high_resolution_clock::now();
    string seq_a = "TCAAGCGTTAGAGAAGTCATTATGTGATAAAAAAATTCAACTTGGTATCAACTTAACTAA";
    string seq_b = "ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCT";

    cout << "Secuencia 01: " << seq_a << endl;
    cout << "Secuencia 02: " << seq_b << endl;

    global_seq(seq_a, seq_b);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    cout << "Tiempo de ejecución: " << duration.count() << " nanosegundos" << endl;
    return 0;
}
