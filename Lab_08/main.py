import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

def couple(x, y):
    pairs = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}
    return -1 if pairs.get(x) == y else 0

def init_matrix(size):
    return np.zeros((size, size), dtype=int)

def fill_matrix(matrix, rna):
    M = len(rna)
    for k in range(1, M):
        for i in range(M - k):
            j = i + k
            left = matrix[i, j - 1]
            down = matrix[i + 1, j]
            diag = matrix[i + 1, j - 1] + couple(rna[i], rna[j])
            rc = np.inf
            for t in range(i + 1, j):
                rc = min(rc, matrix[i, t] + matrix[t + 1, j])
            matrix[i, j] = min(left, down, diag, rc)

def traceback_all(matrix, rna, fold, i, j, results):
    if i < j:
        if i < j - 2 and matrix[i, j] == matrix[i + 1, j - 1] + couple(rna[i], rna[j]):
            fold.append((i, j))
            traceback_all(matrix, rna, fold, i + 1, j - 1, results)
            fold.pop()
        if matrix[i, j] == matrix[i, j - 1]:
            traceback_all(matrix, rna, fold, i, j - 1, results)
        if matrix[i, j] == matrix[i + 1, j]:
            traceback_all(matrix, rna, fold, i + 1, j, results)
        for k in range(i + 1, j):
            if matrix[i, j] == matrix[i, k] + matrix[k + 1, j]:
                traceback_all(matrix, rna, fold, i, k, results)
                traceback_all(matrix, rna, fold, k + 1, j, results)
    else:
        if fold not in results:  # Avoid duplication
            results.append(fold.copy())

def dot_write(rna, fold):
    dot = ['.'] * len(rna)
    for i, j in fold:
        dot[i] = '('
        dot[j] = ')'
    return ''.join(dot)

def draw_rna_structure_circular(rna, fold):
    G = nx.Graph()
    for i in range(len(rna)):
        G.add_node(i, label=rna[i])
    for (i, j) in fold:
        G.add_edge(i, j)

    pos = nx.circular_layout(G)  # Circular layout for the RNA structure
    labels = nx.get_node_attributes(G, 'label')
    nx.draw_networkx_nodes(G, pos, node_color='lightblue', node_size=700)
    nx.draw_networkx_labels(G, pos, labels, font_size=10, font_weight='bold')
    nx.draw_networkx_edges(G, pos, edgelist=fold, edge_color='red', width=2)
    plt.title("Circular RNA Secondary Structure")
    plt.axis('equal')  # Equal aspect ratio
    plt.axis('off')  # Hide the axes
    plt.show()

if __name__ == "__main__":
    rna = "GGAAAUCC"
    matrix = init_matrix(len(rna))
    fill_matrix(matrix, rna)

    results = []
    fold = []
    traceback_all(matrix, rna, fold, 0, len(rna) - 1, results)

    print("RNA SEQUENCE:", rna)
    print("ENERGY MATRIX:")
    for row in matrix:
        print('\t'.join(map(str, row)))

    unique_results = []
    for result in results:
        dot_bracket = dot_write(rna, result)
        if dot_bracket not in unique_results:
            unique_results.append(dot_bracket)
            print("DOT-BRACKET NOTATION:", dot_bracket)
            draw_rna_structure_circular(rna, result)
