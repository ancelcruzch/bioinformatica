#include "clustering.hpp"
#include <iostream>

extern int index_f;
extern int index_c;

void MatrixWithLabels::printMatrix(){
    printMatrix(this->matrix);
}

void MatrixWithLabels::printMatrix(vector<vector<double> >& _matrix){
    // Print column labels
    std::cout << " \t";
    for (const auto& label : labels) {
        std::cout << label << "\t";
    }
    std::cout << std::endl;

    // Print each row with its label
    for (size_t i = 0; i < _matrix.size(); ++i) {
        std::cout << labels[i] << "\t";
        for (const auto& value : _matrix[i]) {
            std::cout << value << "\t";
        }
        std::cout << std::endl;
    }
}

void MatrixWithLabels::makeSymmetrical() {
    makeSymmetrical(this->matrix);
}

void MatrixWithLabels::makeSymmetrical(vector<vector<double> >& _matrix){
    int n = _matrix.size();
    for (int i = 0; i < n; ++i) {
        for (int j = i; j < n; ++j) {
            if(i==j){continue;}
            _matrix[i][j] = _matrix[j][i];
        }
    }
}

double MatrixWithLabels::find_max(){
    double max = DBL_MIN;
    int n = matrix.size();
    for(int i=0; i<n; i++){
        for(int j=i; j<n; j++){
            if(i==j){continue;}
            if(matrix[j][i] > max){
                max = matrix[j][i];
                index_f = j;
                index_c = i;
            }
        }
    }
    return max;
}

double MatrixWithLabels::find_min(){
    double min = DBL_MAX;
    int n = matrix.size();
    for(int i=0; i<n; i++){
        for(int j=i; j<n; j++){
            if(i==j){continue;}
            if(matrix[j][i] < min){
                min = matrix[j][i];
                index_f = j;
                index_c = i;
            }
        }
    }
    return min;
}

void MatrixWithLabels::minimumDistance(int f, int c){
    // Copia de matriz y etiquetas anteriores..
    vector<std::vector<double> > copy_matrix = this->matrix;
    vector<std::string> copy_labels = this->labels;
    makeSymmetrical(copy_matrix);
    cout << "Posicion:  ["<< f << "," << c << "]" << endl;

    int min_aux = min(f, c);
    int max_aux = max(f, c);
    string str_aux = "";

    // Eliminar columnas..
    str_aux += labels[min_aux];
    labels.erase(labels.begin() + min_aux);
    for (size_t i = 0; i < matrix.size(); ++i) {
        matrix[i].erase(matrix[i].begin() + min_aux);
    }
    str_aux += labels[max_aux-1];
    labels.erase(labels.begin() + max_aux-1);
    for (size_t i = 0; i < matrix.size(); ++i) {
        matrix[i].erase(matrix[i].begin() + max_aux-1);
    }
    //Eliminar filas..
    matrix.erase(matrix.begin() + min_aux);
    matrix.erase(matrix.begin() + max_aux-1);

    // Agregar una nueva columna en la posición 'posicion'
    for (int i = 0; i < matrix.size(); ++i) {
        matrix[i].insert(matrix[i].begin() + min_aux, 0.0);
    }

    // Agregar una nueva fila en la posición 'posicion'
    matrix.insert(matrix.begin() + min_aux, std::vector<double>(matrix[0].size(), 0.0));

    // Agregar una nueva etiqueta en la posición 'posicion'
    labels.insert(labels.begin() + min_aux, str_aux);

    int new_n = matrix.size();

    // Actualizar valores nuevos..
    for(int w=0; w<min_aux; w++){
        matrix[min_aux][w] = 0;
        string tmp = labels[w];
        auto it = find(copy_labels.begin(), copy_labels.end(), tmp);
        int index = -1;
        if (it != copy_labels.end()) {
            index = distance(copy_labels.begin(), it);
        }
        matrix[min_aux][w] = min(copy_matrix[min_aux][index], copy_matrix[max_aux][index]);
    }
    for(int w=min_aux+1; w<new_n; w++){
        matrix[w][min_aux] = 0;
        string tmp = labels[w];
        auto it = find(copy_labels.begin(), copy_labels.end(), tmp);
        int index = -1;
        if (it != copy_labels.end()) {
            index = distance(copy_labels.begin(), it);
        }
        matrix[w][min_aux] = min(copy_matrix[min_aux][index], copy_matrix[max_aux][index]);
    }

}

void MatrixWithLabels::maximumDistance(int f, int c){
    // Copia de matriz y etiquetas anteriores..
    vector<std::vector<double> > copy_matrix = this->matrix;
    vector<std::string> copy_labels = this->labels;
    makeSymmetrical(copy_matrix);
    int min_aux = min(f, c);
    int max_aux = max(f, c);
    string str_aux = "";

    // Eliminar columnas..
    str_aux += labels[min_aux];
    labels.erase(labels.begin() + min_aux);
    for (size_t i = 0; i < matrix.size(); ++i) {
        matrix[i].erase(matrix[i].begin() + min_aux);
    }
    str_aux += labels[max_aux-1];
    labels.erase(labels.begin() + max_aux-1);
    for (size_t i = 0; i < matrix.size(); ++i) {
        matrix[i].erase(matrix[i].begin() + max_aux-1);
    }
    //Eliminar filas..
    matrix.erase(matrix.begin() + min_aux);
    matrix.erase(matrix.begin() + max_aux-1);

    // Agregar una nueva columna en la posición 'posicion'
    for (int i = 0; i < matrix.size(); ++i) {
        matrix[i].insert(matrix[i].begin() + min_aux, 0.0);
    }

    // Agregar una nueva fila en la posición 'posicion'
    matrix.insert(matrix.begin() + min_aux, std::vector<double>(matrix[0].size(), 0.0));

    // Agregar una nueva etiqueta en la posición 'posicion'
    labels.insert(labels.begin() + min_aux, str_aux);

    int new_n = matrix.size();

    // Actualizar valores nuevos..
    for(int w=0; w<min_aux; w++){
        matrix[min_aux][w] = 0;
        string tmp = labels[w];
        auto it = find(copy_labels.begin(), copy_labels.end(), tmp);
        int index = -1;
        if (it != copy_labels.end()) {
            index = distance(copy_labels.begin(), it);
        }
        matrix[min_aux][w] = max(copy_matrix[min_aux][index], copy_matrix[max_aux][index]);
    }
    for(int w=min_aux+1; w<new_n; w++){
        matrix[w][min_aux] = 0;
        string tmp = labels[w];
        auto it = find(copy_labels.begin(), copy_labels.end(), tmp);
        int index = -1;
        if (it != copy_labels.end()) {
            index = distance(copy_labels.begin(), it);
        }
        matrix[w][min_aux] = max(copy_matrix[min_aux][index], copy_matrix[max_aux][index]);
    }

}

void MatrixWithLabels::weightedDistance(int f, int c){
    // Copia de matriz y etiquetas anteriores..
    vector<std::vector<double> > copy_matrix = this->matrix;
    vector<std::string> copy_labels = this->labels;
    makeSymmetrical(copy_matrix);
    int min_aux = min(f, c);
    int max_aux = max(f, c);
    string str_aux = "";

    // Eliminar columnas..
    str_aux += labels[min_aux];
    labels.erase(labels.begin() + min_aux);
    for (size_t i = 0; i < matrix.size(); ++i) {
        matrix[i].erase(matrix[i].begin() + min_aux);
    }
    str_aux += labels[max_aux-1];
    labels.erase(labels.begin() + max_aux-1);
    for (size_t i = 0; i < matrix.size(); ++i) {
        matrix[i].erase(matrix[i].begin() + max_aux-1);
    }
    //Eliminar filas..
    matrix.erase(matrix.begin() + min_aux);
    matrix.erase(matrix.begin() + max_aux-1);

    // Agregar una nueva columna en la posición 'posicion'
    for (int i = 0; i < matrix.size(); ++i) {
        matrix[i].insert(matrix[i].begin() + min_aux, 0.0);
    }

    // Agregar una nueva fila en la posición 'posicion'
    matrix.insert(matrix.begin() + min_aux, std::vector<double>(matrix[0].size(), 0.0));

    // Agregar una nueva etiqueta en la posición 'posicion'
    labels.insert(labels.begin() + min_aux, str_aux);

    int new_n = matrix.size();

    // Actualizar valores nuevos..
    for(int w=0; w<min_aux; w++){
        matrix[min_aux][w] = 0;
        string tmp = labels[w];
        auto it = find(copy_labels.begin(), copy_labels.end(), tmp);
        int index = -1;
        if (it != copy_labels.end()) {
            index = distance(copy_labels.begin(), it);
        }
        matrix[min_aux][w] = (copy_matrix[min_aux][index] + copy_matrix[max_aux][index])/2;
    }
    for(int w=min_aux+1; w<new_n; w++){
        matrix[w][min_aux] = 0;
        string tmp = labels[w];
        auto it = find(copy_labels.begin(), copy_labels.end(), tmp);
        int index = -1;
        if (it != copy_labels.end()) {
            index = distance(copy_labels.begin(), it);
        }
        matrix[w][min_aux] = (copy_matrix[min_aux][index] + copy_matrix[max_aux][index])/2;
    }
}
