#include <algorithm>
#include <climits>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>
#include <sstream>
#include <cfloat>
#include <cmath>
#include <numeric>
#include <stdexcept>

#include "clustering.hpp"

int index_f;
int index_c;

int ccc1, ccc2, ccc3;

void readFile(vector<vector<double> > &matrix, vector<string> &labels) {
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
    vector<std::vector<double> > copy_matrix = Mat;
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

    cout << endl;

    do{
        double tmp_min = root.find_min();
        cout << "Nivel K=" << nivel++ << endl;
        cout << "Clusteres mas proximos: (" << root.labels[index_f] << "," << root.labels[index_c] << ") -> " << tmp_min << endl;
        txt_2 << "(" << root.labels[index_f] << "," << root.labels[index_c] << "," << tmp_min << ")" << endl;
        root.setCofeMatrix(root.labels[index_f], root.labels[index_c], tmp_min);
        root.minimumDistance(index_f, index_c);
        cout << "->Matriz quedaria asi: " << endl;
        root.printMatrix();
        cout << "Clusteres: ";
        for(const auto &i: root.labels){
            cout<<i<<", ";
        }
        cout << endl;
        cout << endl;

        /* +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+- */

    }while(root.labels.size()>1);

    cout << "+-+-+-+-+-+-+-+-+-+-+-+ Matriz Cofenetica +-+-+-+-+-+-+-+-+-+-+-+ " << endl;
    root.cofenetica.printMatrix();
    cout << endl;
    double CCC = ccc(copy_matrix, root.cofenetica.matrix);
    ccc1 = CCC;
    cout << endl;
    cout << "CCC: " << CCC << endl;
    return;
}

void distancia_maxima(vector<std::vector<double> > &Mat, vector<string> &Eti, ofstream &txt_2) {
    vector<std::vector<double> > copy_matrix = Mat;
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

    cout << endl;

    do{
        double tmp_min = root.find_min();
        cout << "Nivel K=" << nivel++ << endl;
        cout << "Clusteres mas proximos: (" << root.labels[index_f] << "," << root.labels[index_c] << ") -> " << tmp_min << endl;
        txt_2 << "(" << root.labels[index_f] << "," << root.labels[index_c] << "," << tmp_min << ")" << endl;
        root.setCofeMatrix(root.labels[index_f], root.labels[index_c], tmp_min);
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

    cout << "+-+-+-+-+-+-+-+-+-+-+-+ Matriz Cofenetica +-+-+-+-+-+-+-+-+-+-+-+ " << endl;
    root.cofenetica.printMatrix();
    cout << endl;
    double CCC = ccc(copy_matrix, root.cofenetica.matrix);
    ccc2 = CCC;
    cout << endl;
    cout << "CCC: " << CCC << endl;
    return;
}

void distancia_ponderada(vector<std::vector<double> > &Mat, vector<string> &Eti, ofstream &txt_2) {
    vector<std::vector<double> > copy_matrix = Mat;
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

    cout << endl;

    do{
        double tmp_min = root.find_min();
        cout << "Nivel K=" << nivel++ << endl;
        cout << "Clusteres mas proximos: (" << root.labels[index_f] << "," << root.labels[index_c] << ") -> " << tmp_min << endl;
        txt_2 << "(" << root.labels[index_f] << "," << root.labels[index_c] << "," << tmp_min << ")" << endl;
        root.setCofeMatrix(root.labels[index_f], root.labels[index_c], tmp_min);
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

    cout << "+-+-+-+-+-+-+-+-+-+-+-+ Matriz Cofenetica +-+-+-+-+-+-+-+-+-+-+-+ " << endl;
    root.cofenetica.printMatrix();
    cout << endl;
    double CCC = ccc(copy_matrix, root.cofenetica.matrix);
    ccc3 = CCC;
    cout << endl;
    cout << "CCC: " << CCC << endl;
    return;
}

/* ----------------------------------------------------------------------------------------------------------------- */

void modify_cluster(vector<string>&clusters, string label){
    for(int i=0; i<clusters.size(); i++){
        // Encuentra la posición del substring a eliminar
        size_t pos = clusters[i].find(label);
        if (pos != std::string::npos) {
            clusters[i].erase(pos, label.length());
            clusters.push_back(label);
            return;
        }
    }
}

void metodos_jerarquicos(vector<std::vector<double> > &Mat, vector<string> &Eti) {

    MatrixWithLabels root(Mat, Eti);
    cout << "+-+-+-+-+-+-+-+-+-+-+-+ Metodos Jerarquicos Disociativos +-+-+-+-+-+-+-+-+-+-+-+" << endl;
    root.makeSymmetrical();
    root.printMatrix();

    int n = root.matrix.size();

    // 0. Generamos el cluster inicial..
    string initial_label = "";
    vector<string>clusters;
    for(int i=0; i<n; i++){
        initial_label += string(1, 'A' + i);
    }
    clusters.push_back(initial_label);
    cout << "Clusters: ";
    for(const auto &i : clusters){
        cout << i << endl;
    }

    cout << endl;

    // 1. Calcular distancias de cada individuo..
    cout << "Distancias de cada invidivio al cluster formado.." << endl;
    vector<pair<string, double> >tmp_dist(n, make_pair("", 0.0));
    for(int i=0; i<n; i++){
        double tmp = DBL_MAX;
        for(int j=0; j<n; j++){
            if(i==j){
                continue;
            }
            if(root.matrix[i][j]<tmp){
                tmp = root.matrix[i][j];
            }
        }
        tmp_dist[i].first = root.labels[i];
        tmp_dist[i].second = tmp;
    }
    // Imprimir distancias de cada individuo..
    for(const auto &i:tmp_dist){
        cout << i.first << "  |  " << i.second << endl;
    }
    cout << endl;

    auto max_element_it = max_element(tmp_dist.begin(), tmp_dist.end());
    int max_position = distance(tmp_dist.begin(), max_element_it);
    double max_value = max_element_it->second;
    string max_letter = max_element_it->first;

    /*
    cout << max_value << endl;
    cout << max_letter << endl;
    cout << max_position << endl;
    */

    string new_label = max_letter;

    cout << "New Cluster: " << new_label << endl;
    cout << endl;
    root.dissociativeMethods(tmp_dist, new_label);

    // Modificar cluster | tmp_dist
    cout << clusters[0] << endl;
    modify_cluster(clusters, new_label);
    for(const auto& i:clusters){
        cout << i << ",";
    }
    cout << endl;
    auto new_end = remove_if(tmp_dist.begin(), tmp_dist.end(),[&new_label](const pair<string, double>& p) {
        return p.first == new_label;
    });
    tmp_dist.erase(new_end, tmp_dist.end());
    cout << endl;

    cout << "size clusters: " << clusters.size() << endl;
    cout << "size tmp_dist: " << tmp_dist.size() << endl;

    return;
}

/* ----------------------------------------------------------------------------------------------------------------- */

int main(){
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

    vector<vector<double> >matrix;
    vector<string>labels;
    readFile(matrix, labels);

    // +-+-+ Metodos aglomerativos +-+-+
    //distancia_minima(matrix, labels, dendo_txt);
    //distancia_maxima(matrix, labels, dendo_txt);
    distancia_ponderada(matrix, labels, dendo_txt);

    // +-+-+ Metodo jearquico disociativo +-+-+
    //metodos_jerarquicos(matrix, labels);

    //demo.printMatrix();

    // Restaurar el buffer original de cout
    std::cout.rdbuf(coutbuf);
    outFile.close();

    dendo_txt.close();
    system("/opt/homebrew/anaconda3/bin/python3 main.py");

    return 0;
}
