#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <queue>
#include <mutex>
#include <thread>

using namespace std;

mutex mtx; // Mutex para proteger el acceso a alignments

// Función para imprimir la matriz de alineamiento
void printMatrix(const vector<vector<int>>& matrix) {
    for (const auto& row : matrix) {
        for (int elem : row) {
            cout << elem << " ";
        }
        cout << "\n";
    }
}

// Función para realizar el alineamiento global en paralelo
void global_seq_parallel(const string& seq1, const string& seq2, vector<pair<string, string>>& alignments, int num_threads, ofstream &file) {
    int match = 1;
    int mismatch = -1;
    int gap = -2;
    int n = seq1.size();
    int m = seq2.size();

    // Inicializar la matriz de alineamiento
    vector<vector<int>> matrix(n + 1, vector<int>(m + 1, 0));
    for (int i = 0; i <= n; ++i) {
        matrix[i][0] = i * gap;
    }
    for (int j = 0; j <= m; ++j) {
        matrix[0][j] = j * gap;
    }

    // Llenar la matriz de alineamiento
    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= m; ++j) {
            int diagonal = matrix[i - 1][j - 1] + (seq1[i - 1] == seq2[j - 1] ? match : mismatch);
            int up = matrix[i - 1][j] + gap;
            int left = matrix[i][j - 1] + gap;
            matrix[i][j] = max({diagonal, up, left});
        }
    }

    // Función para hacer el backtracking y encontrar todos los alineamientos óptimos
    auto traceback = [&](int i, int j, string aligned1, string aligned2) {
        queue<tuple<int, int, string, string>> jobs;
        jobs.emplace(i, j, aligned1, aligned2);

        while (!jobs.empty()) {
            tie(i, j, aligned1, aligned2) = jobs.front();
            jobs.pop();

            if (i == 0 && j == 0) {
                reverse(aligned1.begin(), aligned1.end());
                reverse(aligned2.begin(), aligned2.end());
                lock_guard<mutex> lock(mtx);
                alignments.push_back({aligned1, aligned2});
                continue;
            }

            if (i > 0 && matrix[i][j] == matrix[i - 1][j] + gap) {
                jobs.emplace(i - 1, j, aligned1 + seq1[i - 1], aligned2 + '-');
            }
            if (j > 0 && matrix[i][j] == matrix[i][j - 1] + gap) {
                jobs.emplace(i, j - 1, aligned1 + '-', aligned2 + seq2[j - 1]);
            }
            if (i > 0 && j > 0 && matrix[i][j] == matrix[i - 1][j - 1] + (seq1[i - 1] == seq2[j - 1] ? match : mismatch)) {
                jobs.emplace(i - 1, j - 1, aligned1 + seq1[i - 1], aligned2 + seq2[j - 1]);
            }
        }
    };

    // Llamar a la función de traceback en paralelo
    traceback(n, m, "", "");

    cout << "Matriz de alineamiento:\n";
    printMatrix(matrix);

    cout << "Alineamientos óptimos:\n";
    for (const auto& alignment : alignments) {
        cout << alignment.first << "\n";
        cout << alignment.second << "\n\n";
    }
}

int main() {
    auto start = std::chrono::high_resolution_clock::now();
    string seq_a = "TCAAGCGTTATCAAGCGTTATCAAGCGTTATCAAGCGTTATCAAGCGTTATCAAGCGTTATCAAGCGTTATCAAGCGTTA";
    string seq_b = "ATGGAAGCAATCAAGCTCAAGCTCAAGCTCAAGCTCAAGCTCAAGC";

    cout << "Secuencia 01: " << seq_a << endl;
    cout << "Secuencia 02: " << seq_b << endl;

    int num_threads = 4;
    vector<pair<string, string>> alignments;

    global_seq_parallel(seq_a, seq_b, alignments, num_threads);

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    cout << "Tiempo de ejecución: " << duration.count() << " nanosegundos" << endl;

    return 0;
}
