/*
 ESTUDIANTE: ANCEL ALAIN FERNANDO, CRUZ CHAIÑA
 CURSO: BIOINFORMATICA GRUPO-A
 TEMA: LABORATORIO_04
 */

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <cctype>
#include <thread>
#include <mutex>
#include <queue>
#include <chrono>

using namespace std;

string seq_a = "ATTAAAGGTT";
string seq_b = "ATGGAAGCAA";

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

//GLOBAL SECUENCIAL

// Función para hacer el backtracking y encontrar todos los alineamientos óptimos.
void traceback1(vector<vector<int>>& matrix, string& seq1, string& seq2, int i, int j, string& aligned1, string& aligned2, int match, int mismatch, int gap, vector<pair<string, string>>& alignments) {
    if (i == 0 && j == 0) {
        reverse(aligned1.begin(), aligned1.end());
        reverse(aligned2.begin(), aligned2.end());
        alignments.push_back({aligned1, aligned2});
        reverse(aligned1.begin(), aligned1.end());
        reverse(aligned2.begin(), aligned2.end());
        //file << "Alineamiento: " << alignments.size() << endl;
        return;
    }

    if (i > 0 && matrix[i][j] == matrix[i-1][j] + gap) {
        aligned1.push_back(seq1[i-1]);
        aligned2.push_back('-');
        traceback1(matrix, seq1, seq2, i-1, j, aligned1, aligned2, match, mismatch, gap, alignments);
        aligned1.pop_back();
        aligned2.pop_back();
    }
    if (j > 0 && matrix[i][j] == matrix[i][j-1] + gap) {
        aligned1.push_back('-');
        aligned2.push_back(seq2[j-1]);
        traceback1(matrix, seq1, seq2, i, j-1, aligned1, aligned2, match, mismatch, gap, alignments);
        aligned1.pop_back();
        aligned2.pop_back();
    }
    if (i > 0 && j > 0 && matrix[i][j] == matrix[i-1][j-1] + (seq1[i-1] == seq2[j-1] ? match : mismatch)) {
        aligned1.push_back(seq1[i-1]);
        aligned2.push_back(seq2[j-1]);
        traceback1(matrix, seq1, seq2, i-1, j-1, aligned1, aligned2, match, mismatch, gap, alignments);
        aligned1.pop_back();
        aligned2.pop_back();
    }
}

void global_seq1(string seq1, string seq2, ofstream &file) {
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
    file << "Matriz de alineamiento:\n";
    for (auto& row : matrix) {
        for (auto& elem : row) {
            file << elem << " ";
        }
        file << "\n";
    }
    file << "Score: " << matrix[seq_a.size()][seq_b.size()] << endl;
    vector<pair<string, string> > alignments;
    string aligned1, aligned2;
    traceback1(matrix, seq1, seq2, n, m, aligned1, aligned2, match, mismatch, gap, alignments);

    file << "Total alineaminetos: " << alignments.size() << endl;

    file << "Alineamientos óptimos:\n";
    for (const auto& alignment : alignments) {
        file << alignment.first << "\n";
        file << alignment.second << "\n";
        file << "\n";
    }
}

int main_secuencial(){
    ofstream file("output_secuencial.txt");
    auto start = std::chrono::high_resolution_clock::now();
    if(file.is_open()){
        file << "+-+-+-+-+-+-+-+-+-+SECUENCIAS+-+-+-+-+-+-+-+-+-+" << endl;
        file << "Secuencia 01: " << seq_a << endl;
        file << "Secuencia 02: " << seq_b << endl;
        global_seq1(seq_a, seq_b,file);
        file << endl;
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
        file << "Tiempo de ejecución: " << duration.count() << " nanosegundos" << endl;
        file.close();
    }
    return 0;
}


// GLOBAL PARALELO

mutex mtx; // Mutex para proteger el acceso a alignments

// Función para hacer el backtracking y encontrar todos los alineamientos óptimos.
void traceback2(vector<vector<int>>& matrix, string& seq1, string& seq2, int i, int j, string aligned1, string aligned2, int match, int mismatch, int gap, vector<pair<string, string>>& alignments) {
    if (i == 0 && j == 0) {
        reverse(aligned1.begin(), aligned1.end());
        reverse(aligned2.begin(), aligned2.end());
        lock_guard<mutex> lock(mtx);
        alignments.push_back({aligned1, aligned2});
        return;
    }

    vector<thread> threads;

    if (i > 0 && matrix[i][j] == matrix[i-1][j] + gap) {
        threads.emplace_back(traceback2, ref(matrix), ref(seq1), ref(seq2), i-1, j, aligned1 + seq1[i-1], aligned2 + '-', match, mismatch, gap, ref(alignments));
    }
    if (j > 0 && matrix[i][j] == matrix[i][j-1] + gap) {
        threads.emplace_back(traceback2, ref(matrix), ref(seq1), ref(seq2), i, j-1, aligned1 + '-', aligned2 + seq2[j-1], match, mismatch, gap, ref(alignments));
    }
    if (i > 0 && j > 0 && matrix[i][j] == matrix[i-1][j-1] + (seq1[i-1] == seq2[j-1] ? match : mismatch)) {
        threads.emplace_back(traceback2, ref(matrix), ref(seq1), ref(seq2), i-1, j-1, aligned1 + seq1[i-1], aligned2 + seq2[j-1], match, mismatch, gap, ref(alignments));
    }

    for (auto& t : threads) {
        t.join();
    }
}

void global_seq2(string seq1, string seq2, ofstream &file) {
    int match = 1;
    int mismatch = -1;
    int gap = -2;
    int n = seq1.size();
    int m = seq2.size();
    vector<vector<int> > matrix(n + 1, vector<int>(m + 1, 0));
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
    file << "Matriz de alineamiento:\n";
    for (auto& row : matrix) {
        for (auto& elem : row) {
            file << elem << " ";
        }
        file << "\n";
    }

    vector<pair<string, string>> alignments;
    string aligned1, aligned2;
    traceback2(matrix, seq1, seq2, n, m, aligned1, aligned2, match, mismatch, gap, alignments);

    file << "Alineamientos óptimos:\n";
    for (const auto& alignment : alignments) {
        file << alignment.first << "\n";
        file << alignment.second << "\n";
        file << "\n";
    }
}

// Función para realizar el alineamiento global en paralelo
void global_seq3(const string& seq1, const string& seq2, vector<pair<string, string>>& alignments, int num_threads, ofstream &file) {
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

    file << "Matriz de alineamiento:\n";
    for (const auto& row : matrix) {
        for (int elem : row) {
            file << elem << " ";
        }
        file << "\n";
    }
    file << "Total alineaminetos: " << alignments.size() << endl;
    file << "Score: " << matrix[seq_a.size()][seq_b.size()] << endl;
    file << "Alineamientos óptimos:\n";

    for (const auto& alignment : alignments) {
        file << alignment.first << "\n";
        file << alignment.second << "\n\n";
    }
}


int main_paralelo(){
    ofstream file("output_paralelo.txt");
    auto start = std::chrono::high_resolution_clock::now();
    if(file.is_open()){
        file << "+-+-+-+-+-+-+-+-+-+SECUENCIAS+-+-+-+-+-+-+-+-+-+" << endl;
        file << "Secuencia 01: " << seq_a << endl;
        file << "Secuencia 02: " << seq_b << endl;
        global_seq2(seq_a, seq_b,file);
        file << endl;
        file.close();
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    cout << "Tiempo de ejecución: " << duration.count() << " nanosegundos" << endl;
    return 0;
}

int main_paralelo2(){
    vector<pair<string, string>> alignments;
    ofstream file("output_paralelo2.txt");
    auto start = std::chrono::high_resolution_clock::now();
    if(file.is_open()){
        file << "+-+-+-+-+-+-+-+-+-+SECUENCIAS+-+-+-+-+-+-+-+-+-+" << endl;
        file << "Secuencia 01: " << seq_a << endl;
        file << "Secuencia 02: " << seq_b << endl;
        global_seq3(seq_a, seq_b,alignments,4,file);
        file << endl;
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
        file << "Tiempo de ejecución: " << duration.count() << " nanosegundos" << endl;
        file.close();
    }
    return 0;
}

int main(){
    main_secuencial();
    //main_paralelo();
    main_paralelo2();
}
