#include <iostream>
#include <vector>
#include <string>

bool isAllGaps(const std::vector<std::string>& alignments, size_t index) {
    for (const auto& alignment : alignments) {
        if (index >= alignment.size() || alignment[index] != '-') {
            return false;
        }
    }
    return true;
}

int main() {
    // Inicializa el vector con los alineamientos
    std::vector<std::string> multiple_alignment = {
        "A-TTGCCATT",
        "A-TGGCCATT",
        "A-TC-CAATTTT",
        "A-TCTTC-TT",
        "A-CTGACC--"
    };

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

    // Imprime el alineamiento múltiple
    for (const auto& alignment : multiple_alignment) {
        std::cout << alignment << std::endl;
    }

    return 0;
}
