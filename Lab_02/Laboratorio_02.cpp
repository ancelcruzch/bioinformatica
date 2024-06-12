/*
 ESTUDIANTE: ANCEL ALAIN FERNANDO, CRUZ CHAIÑA
 CURSO: BIOINFORMATICA GRUPO-A
 TEMA: LABORATORIO_02
 */

#include <iostream>
#include <pthread.h>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <cctype>

using namespace std;

// Pregunta 01
// El primer parametro sera el substring a tomar en cuenta, mientras que el segundo sera el string mayor.
bool es_substring(string sub_str, string str){
    return str.find(sub_str) != string::npos;
}

// Pregunta 02
// Funcion de score teniendo en cuenta, si son iguales valor de +1 y si son diferentes -2.
int def_score(string str_1, string str_2)
{
    int score = 0;
    int min_size = min(str_1.size(), str_2.size());
    int max_size = max(str_1.size(), str_2.size());
    for (int i = 0; i < min_size; i++){
        if (str_1[i] == str_2[i]){
            score += 1;
        }else{
            score -= 2;
        }
    }
    // cout << score << endl;
    score -= (2) * (max_size - min_size); // -2 para el resto de las secuencias..
    cout << "Score: " << score << endl;
    return score;
}

int main(){
    return 0;
}
