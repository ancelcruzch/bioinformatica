#ifndef CLUSTERING_HPP
#define CLUSTERING_HPP

#include <algorithm>
#include <climits>
#include <string>
#include <type_traits>
#include <vector>
#include <cfloat>

using namespace std;

// ESTRUCTURA DE LA MATRIZ..

struct MatrixWithLabels {
    vector<std::vector<double> > matrix;
    vector<std::string> labels;
    // Constructor..
    MatrixWithLabels(const std::vector<std::vector<double> >& mat, const std::vector<std::string>& lab) : matrix(mat), labels(lab) {}
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
    void readFile();
};


#endif
