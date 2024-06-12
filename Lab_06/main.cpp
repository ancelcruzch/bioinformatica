#include <algorithm>
#include <climits>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <type_traits>
#include <vector>
#include <sstream>
#include <cfloat>

#include "clustering.hpp"

using namespace std;

int index_f, index_c;

void readFile(vector<vector<double> >& matrix, vector<string>& labels) {
    ifstream file("input.txt");
    string line;
    int rowCount = 0;
    int maxCols = 0;
    while (getline(file, line)) {
        vector<double> row;
        istringstream ss(line);
        double value;
        while (ss >> value) {
            row.push_back(value);
        }
        matrix.push_back(row);
        rowCount++;
        // Actualizar el número máximo de columnas si es necesario
        if (row.size() > maxCols) {
            maxCols = row.size();
        }
    }

    // Asegurarse de que la matriz sea cuadrada agregando ceros si es necesario
    for (auto& row : matrix) {
        while (row.size() < maxCols) {
            row.push_back(0.0); // Agregar ceros hasta que tenga el mismo número de columnas que maxCols
        }
    }

    // Generar etiquetas para cada fila
    for (int i = 0; i < rowCount; ++i) {
        // string indexLabel = "A" + to_string(i + 1);
        string indexLabel(1, 'A' + i);
        labels.push_back(indexLabel);
    }


    file.close();
}

void distancia_minima(vector<std::vector<double> > &Mat, vector<string> &Eti, ofstream &txt_2) {

    MatrixWithLabels root(Mat, Eti);
    cout << "+-+-+-+-+-+-+-+-+-+-+-+ Distancia Minima +-+-+-+-+-+-+-+-+-+-+-+" << endl;
    root.printMatrix();
    cout << "Maximo: " << root.find_max() << endl;
    cout << "Minimo: " << root.find_min() << endl;
    cout << "Clusteres: ";
    for(const auto &i: root.labels){
        cout<<i<<", ";
        txt_2<<i<<", ";
    }
    cout << endl;
    txt_2 << endl;

    // ITERACION..
    int nivel = 1;
    int count = 5;

    cout << endl;

    do{
        double tmp_min = root.find_min();
        cout << "Nivel K=" << nivel++ << endl;
        cout << "Clusteres mas proximos: (" << root.labels[index_f] << "," << root.labels[index_c] << ") -> " << tmp_min << endl;
        txt_2 << "(" << root.labels[index_f] << "," << root.labels[index_c] << "," << tmp_min << ")" << endl;
        root.minimumDistance(index_f, index_c);
        cout << "->Matriz quedaria asi: " << endl;
        root.printMatrix();
        cout << "Clusteres: ";
        for(const auto &i: root.labels){
            cout<<i<<", ";
        }
        cout << endl;
        cout << endl;
    }while(root.labels.size()>1);

    return;
}

void distancia_maxima(vector<std::vector<double> > &Mat, vector<string> &Eti, ofstream &txt_2) {

    MatrixWithLabels root(Mat, Eti);
    cout << "+-+-+-+-+-+-+-+-+-+-+-+ Distancia Maxima +-+-+-+-+-+-+-+-+-+-+-+" << endl;
    root.printMatrix();
    cout << "Maximo: " << root.find_max() << endl;
    cout << "Minimo: " << root.find_min() << endl;
    cout << "Clusteres: ";
    for(const auto &i: root.labels){
        cout<<i<<", ";
        txt_2<<i<<", ";
    }
    cout << endl;
    txt_2 << endl;

    // ITERACION..
    int nivel = 1;
    int count = 5;

    cout << endl;

    do{
        double tmp_min = root.find_min();
        cout << "Nivel K=" << nivel++ << endl;
        cout << "Clusteres mas proximos: (" << root.labels[index_f] << "," << root.labels[index_c] << ") -> " << tmp_min << endl;
        txt_2 << "(" << root.labels[index_f] << "," << root.labels[index_c] << "," << tmp_min << ")" << endl;
        root.maximumDistance(index_f, index_c);
        cout << "->Matriz quedaria asi: " << endl;
        root.printMatrix();
        cout << "Clusteres: ";
        for(const auto &i: root.labels){
            cout<<i<<", ";
        }
        cout << endl;
        cout << endl;
    }while(root.labels.size()>1);

    return;
}

void distancia_ponderada(vector<std::vector<double> > &Mat, vector<string> &Eti, ofstream &txt_2) {

    MatrixWithLabels root(Mat, Eti);
    cout << "+-+-+-+-+-+-+-+-+-+-+-+ Distancia Ponderada +-+-+-+-+-+-+-+-+-+-+-+" << endl;
    root.printMatrix();
    cout << "Maximo: " << root.find_max() << endl;
    cout << "Minimo: " << root.find_min() << endl;
    cout << "Clusteres: ";
    for(const auto &i: root.labels){
        cout<<i<<", ";
        txt_2<<i<<", ";
    }
    cout << endl;
    txt_2 << endl;

    // ITERACION..
    int nivel = 1;
    int count = 5;

    cout << endl;

    do{
        double tmp_min = root.find_min();
        cout << "Nivel K=" << nivel++ << endl;
        cout << "Clusteres mas proximos: (" << root.labels[index_f] << "," << root.labels[index_c] << ") -> " << tmp_min << endl;
        txt_2 << "(" << root.labels[index_f] << "," << root.labels[index_c] << "," << tmp_min << ")" << endl;
        root.weightedDistance(index_f, index_c);
        cout << "->Matriz quedaria asi: " << endl;
        root.printMatrix();
        cout << "Clusteres: ";
        for(const auto &i: root.labels){
            cout<<i<<", ";
        }
        cout << endl;
        cout << endl;
    }while(root.labels.size()>1);

    return;
}



int main() {

    // Crear y abrir un archivo de texto
    ofstream outFile("clustering.txt");

    // Verificar si el archivo se abrió correctamente
    if (!outFile) {
        cerr << "No se pudo abrir el archivo para escribir." << std::endl;
        return 1;
    }

    // Redirigir cout a outFile
    std::streambuf* coutbuf = std::cout.rdbuf(); // Guardar el buffer original de cout
    std::cout.rdbuf(outFile.rdbuf()); // Redirigir cout a outFile


    ofstream dendo_txt("dendo.txt");

    vector<string>labels = {"A","B","C","D","E","F","G"};
    vector<vector<double> > matriz = {
        {0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00},
        {2.15, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00},
        {0.70, 1.53, 0.00, 0.00, 0.00, 0.00, 0.00},
        {1.07, 1.14, 0.43, 0.00, 0.00, 0.00, 0.00},
        {0.85, 1.38, 0.21, 0.29, 0.00, 0.00, 0.00},
        {1.16, 1.01, 0.55, 0.22, 0.41, 0.00, 0.00},
        {1.56, 2.83, 1.86, 2.04, 2.02, 2.05, 0.00}
    };

    vector<vector<double> > matriz2;
    vector<string> labels2;
    readFile(matriz2, labels2);

    //distancia_minima(matriz2, labels2,dendo_txt);
    //distancia_maxima(matriz2, labels2,dendo_txt);
    distancia_ponderada(matriz2, labels2, dendo_txt);


    cout << "Longitud: " << matriz2.size() << endl;

    // Restaurar el buffer original de cout
    std::cout.rdbuf(coutbuf);
    outFile.close();

    dendo_txt.close();

    system("python3 main.py");


    return 0;
}
