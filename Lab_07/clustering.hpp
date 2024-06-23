#ifndef CLUSTERING_HPP
#define CLUSTERING_HPP

#include <algorithm>
#include <climits>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>
#include <cfloat>

using namespace std;

struct Cofenetica{
    vector<vector<double> > matrix;
    vector<string> labels;
    // Constructor que copia matrix y labels de MatrixWithLabels..
    Cofenetica(const vector<vector<double> > mat, const vector<string> lab) : matrix(mat), labels(lab) {}
    // Funciones..
    void printMatrix();
    int indexLabel(string);
    void updateValues(string,string,double);
};

// ESTRUCTURA DE LA MATRIZ..
struct MatrixWithLabels {
    vector<std::vector<double> > matrix;
    vector<std::string> labels;
    Cofenetica cofenetica;
    // Constructor..
    MatrixWithLabels(const vector<vector<double> >&mat, const vector<string> &lab);
    // Funciones
    void printMatrix();
    void printMatrix(vector<vector<double> >& _matrix);
    void makeSymmetrical();
    void makeSymmetrical(vector<vector<double> >& _matrix);
    double find_max();
    double find_min();
    void minimumDistance(int f, int c);
    void maximumDistance(int f, int c);
    void weightedDistance(int f, int c);
    double main_distance(string,string);
    void setCofeMatrix(string label_x, string label_y, double value);
    void dissociativeMethods(vector<pair<string, double> >&, string &);
};

double mean(const vector<double> &v);
double variance(const vector<double> &v);
double covariance(const vector<double> &v1, const vector<double> &v2);
double ccc(const vector<vector<double> > &dist_matrix, const vector<vector<double> > &cof_matrix);

#endif
