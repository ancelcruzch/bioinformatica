#include "clustering.hpp"
#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include <numeric>

extern int index_f;
extern int index_c;

MatrixWithLabels::MatrixWithLabels(const vector<vector<double> >& mat, const vector<string>& lab):
    matrix(mat),
    labels(lab),
    cofenetica(mat, lab)
{}

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

double MatrixWithLabels::main_distance(string label1, string label2){
    int n = this->matrix.size();
    auto it1 = find(this->labels.begin(), this->labels.end(), label1);
    auto it2 = find(this->labels.begin(), this->labels.end(), label2);
    if (it1 == this->labels.end() || it2 == this->labels.end()) {
        throw invalid_argument("Etiqueta no encontrada");
    }
    // Calcula los índices
    size_t index1 = distance(this->labels.begin(), it1);
    size_t index2 = distance(this->labels.begin(), it2);
    return this->matrix[index1][index2];
}

void MatrixWithLabels::setCofeMatrix(string label_x, string label_y, double value){
    this->cofenetica.updateValues(label_x, label_y, value);
}

void MatrixWithLabels::dissociativeMethods(vector<pair<string, double> > &main_dist, string &label){
    struct Individual{
        string indiv;
        double distanciaPrinc;
        double distNuevo;
        double difer;
        Individual(string i, double dp, double dn, double df) : indiv(i), distanciaPrinc(dp), distNuevo(dn), difer(df) {}
    };
    vector<Individual>distance_matrix;

    cout << "Invid" << "\t" << "DistP" << "\t" << "DistN" << "\t" << "Difer." << endl;
    for(int i=0;i<main_dist.size(); i++){
        if(main_dist[i].first == label){ // obviar el label ya ingresado
            continue;
        }
        //cout << main_dist[i].first << "\t" << main_dist[i].second << "\t" << main_distance(main_dist[i].first, label) << endl;
        double tmp_main_dist = main_distance(main_dist[i].first, label);
        double tmp_dif = main_dist[i].second - tmp_main_dist;
        distance_matrix.push_back(Individual(main_dist[i].first, main_dist[i].second, tmp_main_dist, tmp_dif));
        cout << main_dist[i].first  << "\t" << main_dist[i].second << "\t" << tmp_main_dist << "\t" << tmp_dif << endl;
    }

    return;
}

/* ----------------------------------------------------------------------------------------------------------------- */

void Cofenetica::printMatrix(){
    // Print column labels
    std::cout << " \t";
    for (const auto& label : labels) {
        std::cout << label << "\t";
    }
    std::cout << std::endl;

    // Print each row with its label
    for (size_t i = 0; i < matrix.size(); ++i) {
        std::cout << labels[i] << "\t";
        for (const auto& value : matrix[i]) {
            std::cout << value << "\t";
        }
        std::cout << std::endl;
    }
}


int Cofenetica::indexLabel(string label){
    for(int i=0; i<this->labels.size(); i++){
        if(label == this->labels[i]){
            return i;
        }
    }
    return -1; // en caso no se encuentre..
}

void Cofenetica::updateValues(string label_x, string label_y, double value){
    for(int i=0; i<label_x.size(); i++){
        for(int j=0; j<label_y.size(); j++){
            string tmp_x(1, label_x[i]);
            string tmp_y(1, label_y[j]);
            int index_x = indexLabel(tmp_x);
            int index_y = indexLabel(tmp_y);
            this->matrix[index_x][index_y] = value;
            this->matrix[index_y][index_x] = value;
        }
    }
}

/* ----------------------------------------------------------------------------------------------------------------- */

double mean(const vector<double>& v) {
    return accumulate(v.begin(), v.end(), 0.0) / v.size();
}

double variance(const vector<double>& v) {
    double m = mean(v);
    double var = 0.0;
    for (double value : v) {
        var += (value - m) * (value - m);
    }
    return var / (v.size() - 1);
}

double covariance(const vector<double>& v1, const vector<double>& v2) {
    double mean1 = mean(v1);
    double mean2 = mean(v2);
    double cov = 0.0;
    for (size_t i = 0; i < v1.size(); ++i) {
        cov += (v1[i] - mean1) * (v2[i] - mean2);
    }
    return cov / (v1.size() - 1);
}

// Función para calcular el coeficiente de correlación cofenética (CCC)
double ccc(const vector<vector<double> > &dist_matrix, const vector<vector<double> > &cof_matrix) {
    // Verificar que las matrices no estén vacías y tengan el mismo tamaño
    if (dist_matrix.empty() || cof_matrix.empty() || dist_matrix.size() != cof_matrix.size()) {
        throw invalid_argument("Las matrices deben tener el mismo tamaño y no deben estar vacías.");
    }

    // Paso 2: Calcular las distancias y las cofenéticas como vectores
    vector<double> dist_vector;
    vector<double> cof_vector;

    for (size_t i = 0; i < dist_matrix.size(); ++i) {
        for (size_t j = 0; j < dist_matrix[i].size(); ++j) {
            if (dist_matrix[i][j] != 0) { // Filtrar los ceros
                dist_vector.push_back(dist_matrix[i][j]);
                cof_vector.push_back(cof_matrix[i][j]);
            }
        }
    }

    // Depuración: Imprimir los vectores de distancias y cofenéticas
    cout << "Distancias: ";
    for (double d : dist_vector) {
        cout << d << " ";
    }
    cout << endl;

    cout << "Cofenéticas: ";
    for (double c : cof_vector) {
        cout << c << " ";
    }
    cout << endl;

    // Paso 3: Calcular la covarianza y las varianzas necesarias
    double cov_dist_cof = covariance(dist_vector, cof_vector);
    double var_dist = variance(dist_vector);
    double var_cof = variance(cof_vector);

    // Depuración: Imprimir las varianzas y covarianza
    cout << "Covarianza: " << cov_dist_cof << endl;
    cout << "Varianza de distancias: " << var_dist << endl;
    cout << "Varianza de cofenéticas: " << var_cof << endl;

    // Verificación de varianzas cero para evitar divisiones por cero
    if (var_dist == 0 || var_cof == 0) {
        cerr << "Error: Varianza cero, no se puede calcular CCC." << endl;
        return nan(""); // Devuelve nan si alguna varianza es cero
    }

    // Paso 4: Calcular el coeficiente de correlación cofenética (CCC)
    double CCC = cov_dist_cof / sqrt(var_dist * var_cof);

    return CCC;
}
