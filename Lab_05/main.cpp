/*
 ESTUDIANTE: ANCEL ALAIN FERNANDO, CRUZ CHAIÑA
 CURSO: BIOINFORMATICA GRUPO-A
 TEMA: LABORATORIO_05
*/

#include <climits>
#include <cstddef>
#include <iostream>
#include <type_traits>
#include <utility>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <cctype>
#include <string>

using namespace std;

/* +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+- */

// Variables Globales..
vector<string>multiple_alignment; // Heuristica para alineamiento de multiples secuencias..
vector<string>gaps;

/* +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+- */

// Funciones adicionales de depuracion..
void printMatrix(vector<vector<int> > &matrix) {
    for (auto& row : matrix) {
        for (auto& elem : row) {
            cout << elem << " ";
        }
        cout << "\n";
    }
}

// Función para hacer el backtracking y encontrar todos los alineamientos óptimos..
void traceback(vector<vector<int> >& matrix, string& seq1, string& seq2, int i, int j, string& aligned1, string& aligned2, int match, int mismatch, int gap, vector<pair<string, string>>& alignments) {
    if (i == 0 && j == 0) {
        reverse(aligned1.begin(), aligned1.end());
        reverse(aligned2.begin(), aligned2.end());
        alignments.push_back({aligned1, aligned2});
        reverse(aligned1.begin(), aligned1.end());
        reverse(aligned2.begin(), aligned2.end());
        return;
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

// Algoritmo de Needleman-Wunsch (alineamiento global)..
int needleman_wunsch(string seq1, string seq2, bool show_matrix=true, bool show_alignments=true) {
    int score = 0;
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

    if(show_matrix){
        cout << "-> Secuencias: " << endl;
        cout << seq1 << endl;
        cout << seq2 << endl;
        cout << "Matriz de alineamiento:\n";
        printMatrix(matrix);
    }

    if(show_alignments){
        vector<pair<string, string> > alignments;
        string aligned1, aligned2;
        traceback(matrix, seq1, seq2, n, m, aligned1, aligned2, match, mismatch, gap, alignments);
        cout << "Alineamientos óptimos: " << alignments.size() << endl;
        for (const auto& alignment : alignments) {
            cout << alignment.first << "\n";
            cout << alignment.second << "\n";
            cout << "\n";
        }
        // Incluyendo secuencia optimas en la vector de alineamientos (global)
        multiple_alignment.push_back(alignments[0].second);
        // incluyendo vectores con gaps..
        gaps.push_back(alignments[0].first);

    }

    return matrix[seq1.size()][seq2.size()];
}

/* +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+- */




/* +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+- */

vector<vector<int> > star_alignment(vector<string> secuencias, int &estrella){
    int n = secuencias.size();
    vector<vector<int> > matriz_scores(n, vector<int>(n, 0));
    vector<string>multiple_alignment;
    for(int i=0; i<n; i++){
        for(int j=i; j<n; j++){
            if(i==j){
                matriz_scores[i][j] = 0;
                continue;
            }
            int score_tmp = needleman_wunsch(secuencias[i], secuencias[j], false, false);
            matriz_scores[i][j] = score_tmp;
            matriz_scores[j][i] = score_tmp;
        }
    }

    cout << "Matriz de scores por pares: " << endl;
    printMatrix(matriz_scores);

    // Calculando estrella..
    vector<int>sum_row_tmp(n,0);
    for(int i=0; i<n; i++){
        int aux = 0;
        for(int j=0; j<n; j++){
            aux += matriz_scores[i][j];
        }
        sum_row_tmp[i] = aux;
        cout << "S" << i+1 << ": "  << sum_row_tmp[i] << endl;
    }

    auto maxElementIter = max_element(sum_row_tmp.begin(), sum_row_tmp.end());
    int maxIndex = distance(sum_row_tmp.begin(), maxElementIter);

    estrella = maxIndex;
    return matriz_scores;
}


/* +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+- */

bool isAllGaps(const std::vector<std::string>& alignments, size_t index) {
    for (const auto& alignment : alignments) {
        if (index >= alignment.size() || alignment[index] != '-') {
            return false;
        }
    }
    return true;
}

void function_star_alignment(){
    // Secuencias de ejemplo (diapositivas)
    vector<string> secuencias = {
        "TCAAGCGTTA",
        "GAGAAGTCAT",
        "TATGTGATAA",
        "AAAAATTCAA",
        "CTTGGTATCA"
    };

    cout << "+-+-+-+-+-+-+-STAR-ALIGNMENT-+-+-+-+-+-+-+" << endl;
    for(const auto& i : secuencias){
        cout << i << endl;
    }
    cout << endl;

    int estrella = INT_MAX;
    star_alignment(secuencias, estrella);
    cout << "Estrella: " << secuencias[estrella] << endl;
    cout << endl;

    // Agregamos la estrella al vector de alineamientos..
    multiple_alignment.push_back(secuencias[estrella]);

    for(int i=0; i<secuencias.size(); i++){
        if(i==estrella){
            continue;
        }
        needleman_wunsch(secuencias[estrella], secuencias[i]);
    }

    cout << "CONSTRUYENDO ALINEAMIENTO MULTIPLE: " << endl;
    for(const auto& i : multiple_alignment){
        cout << i << endl;
    }

    cout << endl;

    for(int i=0; i<gaps.size(); i++){
        cout << "gaps: " << gaps[i] << endl;
        vector<int> indices;

        for (int j = 0; j < gaps[i].size(); ++j) {
            if (gaps[j] == "-") {
                indices.push_back(j);
            }
        }
        cout << "Indices: " << endl;
        for (int indice : indices) {
            cout << indice << " ";
        }

    }

    // Encuentra la longitud del alineamiento más largo
    size_t max_length = 0;
    for (const auto& alignment : multiple_alignment) {
        if (alignment.size() > max_length) {
            max_length = alignment.size();
        }
    }

    // Añade guiones al final de los alineamientos más cortos
    for (auto& alignment : multiple_alignment) {
        while (alignment.size() < max_length) {
            alignment += '-';
        }
    }

    // Elimina las columnas de solo guiones
    for (size_t i = 0; i < max_length; ++i) {
        if (isAllGaps(multiple_alignment, i)) {
            for (auto& alignment : multiple_alignment) {
                alignment.erase(i, 1);
            }
            --i;
            --max_length;
        }
    }

    cout << endl;
    cout << "Alineamiento en Estrella: " << endl;

    // Imprime el alineamiento..
    for (const auto& alignment : multiple_alignment) {
        std::cout << alignment << std::endl;
    }

}

/* +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+- */

int main(){

    // Crear y abrir un archivo de texto
    ofstream outFile("star_alignment.txt");

    // Verificar si el archivo se abrió correctamente
    if (!outFile) {
        cerr << "No se pudo abrir el archivo para escribir." << std::endl;
        return 1;
    }

    // Redirigir cout a outFile
    std::streambuf* coutbuf = std::cout.rdbuf(); // Guardar el buffer original de cout
    std::cout.rdbuf(outFile.rdbuf()); // Redirigir cout a outFile

    function_star_alignment();

    // Restaurar el buffer original de cout
    std::cout.rdbuf(coutbuf);
    outFile.close();

    return 0;
}
