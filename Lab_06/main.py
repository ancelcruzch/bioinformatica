import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram
import numpy as np

elementos_originales = []

def leer_uniones_desde_txt(ruta_archivo):
    global elementos_originales  # Para modificar la variable global

    with open(ruta_archivo, 'r') as archivo:
        lineas = archivo.readlines()

    # Leer y procesar la primera línea para obtener los elementos originales
    primera_linea = lineas[0].strip().replace(' ', '').replace(',', '')
    elementos_originales = list(primera_linea)

    # Ignorar la primera línea y procesar las demás líneas
    lineas = lineas[1:]

    uniones = []
    for linea in lineas:
        # Remover paréntesis y espacios, luego dividir por coma
        linea = linea.strip().replace('(', '').replace(')', '')
        partes = linea.split(',')

        # Asegurarse de que hay tres partes
        if len(partes) == 3:
            elem1, elem2, peso = partes[0].strip(), partes[1].strip(), float(partes[2].strip())
            uniones.append((elem1, elem2, peso))

    return uniones

# Ruta del archivo .txt
ruta_archivo = 'dendo.txt'
uniones = leer_uniones_desde_txt(ruta_archivo)
print(elementos_originales)
print(uniones)

# Definir las uniones y pesos
uniones = leer_uniones_desde_txt(ruta_archivo)

# Crear un diccionario para rastrear los índices de los elementos
indice_elementos = {elemento: i for i, elemento in enumerate(elementos_originales)}
nuevo_indice = len(elementos_originales)

# Inicializar las listas para las uniones y los pesos
linkage_matrix = []

# Procesar cada unión
for elem1, elem2, peso in uniones:
    if elem1 in indice_elementos:
        idx1 = indice_elementos[elem1]
    else:
        idx1 = indice_elementos[elem1] = nuevo_indice
        nuevo_indice += 1

    if elem2 in indice_elementos:
        idx2 = indice_elementos[elem2]
    else:
        idx2 = indice_elementos[elem2] = nuevo_indice
        nuevo_indice += 1

    linkage_matrix.append([idx1, idx2, peso, len(elem1) + len(elem2)])

# Convertir a numpy array para usar con dendrogram
linkage_matrix = np.array(linkage_matrix)

# Crear el dendrograma
plt.figure(figsize=(10, 7))
dendrogram(linkage_matrix, labels=elementos_originales)
plt.title("Dendrograma")
plt.xlabel("Elementos")
plt.ylabel("Distancia")
plt.show()
