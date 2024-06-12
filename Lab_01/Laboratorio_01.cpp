/*
 ESTUDIANTE: ANCEL ALAIN FERNANDO, CRUZ CHAIÑA
 CURSO: BIOINFORMATICA GRUPO-A
 TEMA: LABORATORIO_01
 */

#include <iostream>
#include <vector>
#include <string>
#include <fstream>

using namespace std;

//Leer cadenas desde un txt..
vector<string> leerCadenas(const string& nombreArchivo) {
    vector<string> cadenas;
    ifstream archivo(nombreArchivo);
    string linea;
    if (archivo.is_open()) {
        while (getline(archivo, linea)) {
            cadenas.push_back(move(linea));
        }
        archivo.close();
    } else {
        cout << "No se pudo abrir el archivo." << endl;
    }

    return cadenas;
}

//Funcion para identificar si es ADN ARN o Proteina..
string identificar_tipo_cadena(const string& cadena) {
    bool es_adn = true;
    for (char base : cadena) {
        if (toupper(base) != 'A' && toupper(base) != 'C' && toupper(base) != 'G' && toupper(base) != 'T') {
            es_adn = false;
            break;
        }
    }
    if (es_adn) {
        return "ADN";
    }
    bool es_arn = true;
    for (char base : cadena) {
        if (toupper(base) != 'A' && toupper(base) != 'C' && toupper(base) != 'G' && toupper(base) != 'U') {
            es_arn = false;
            break;
        }
    }
    if (es_arn) {
        return "ARN";
    }
    // Si no es ADN ni ARN, asumimos que es una cadena de proteínas
    return "Proteína";
}

//Funcion para la transcripcion de cadenas de ADN a ARN..
string transcripcion_adn_a_arn(const string& adn) {
    string arn = adn;
    for (char& base : arn) {
        if (base == 'T') {
            base = 'U';
        }
    }
    return arn;
}

int main() {
    vector<string>cadenas = { //cadenas conjuntas
        "ATGCGATACGCTTGA",
        "CCGTAAGTTCAACGG",
        "AUGCGAUACGCUUGA",
        "CCGUUAGUUCUACGG",
        "Met-Arg-Tyr-Ala-Stop",
        "Pro-Lys-Phe-Gln-Arg"
    };
    
    vector<string>cadenas_adn = {
        "ATGCGATACGCTTGACCTGATAGCTACGATCGATCGATACGATCGATA",
        "CCGTAAGTTCAACGGCTGATCGTAGCTAGCTAGCTAGCTAGCTAGCTA",
        "TTGACGTATGCTACGATCGATACGATCGATACGATCGATACGATCGAT",
        "GCTATCGATCGATACGATCGATACGATCGATACGATCGATACGATCGA",
        "CGTACGATCGATACGATCGATACGATCGATACGATCGATACGATCGAT",
        "TACGATCGATACGATCGATACGATCGATACGATCGATACGATCGATAC",
        "ATCGATACGATCGATACGATCGATACGATCGATACGATCGATACGATC",
        "CGATACGATCGATACGATCGATACGATCGATACGATCGATACGATCGA",
        "GATACGATCGATACGATCGATACGATCGATACGATCGATACGATCGAT",
        "TACGATCGATACGATCGATACGATCGATACGATCGATACGATCGATAC"
    };
    
    // IMPRIMIR INFORMACION EN CONSOLA.
    cout << "+-+-+-+-+-+-+-+-+-+IDENTIFICAR TIPO DE CADENA+-+-+-+-+-+-+-+-+-+" << endl;
    for(int i=0; i<cadenas.size(); i++){
        cout << i+1 << ". " << cadenas[i] << " : " << identificar_tipo_cadena(cadenas[i]) << endl;
    }
    
    
    cout << "+-+-+-+-+-+-+-+-+-+TRANSCRIPCION DE ADN A ARN+-+-+-+-+-+-+-+-+-+" << endl;
    for(int i=0; i<cadenas_adn.size(); i++){
        cout << i+1 << ". ADN: " << cadenas_adn[i] << " -> ARN: " << transcripcion_adn_a_arn(cadenas_adn[i]) << endl;
    }
    
    //GUARDAR INFORMACION EN UN ARCHIVO TXT.
    string nombre_archivo = "laboratorio_01.txt";
    ofstream archivo(nombre_archivo);
    if (archivo.is_open()) {
        archivo << "+-+-+-+-+-+-+-+-+-+IDENTIFICAR TIPO DE CADENA+-+-+-+-+-+-+-+-+-+" << endl;
        for(int i=0; i<cadenas.size(); i++){
            archivo << i+1 << ". " << cadenas[i] << " : " << identificar_tipo_cadena(cadenas[i]) << endl;
        }
        archivo << endl;
        archivo << "+-+-+-+-+-+-+-+-+-+TRANSCRIPCION DE ADN A ARN+-+-+-+-+-+-+-+-+-+" << endl;
        for(int i=0; i<cadenas_adn.size(); i++){
            archivo << "ADN: " << cadenas_adn[i] << endl;
            archivo << "ARN: " << transcripcion_adn_a_arn(cadenas_adn[i]) << endl;
            archivo << endl;
        }
    }
    
    return 0;
}
