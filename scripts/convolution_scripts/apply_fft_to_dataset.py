import pandas as pd
import sys
import numpy as np

#Nota, todas deben tener el mismo numero de columnas
#NOTE: el 30 es solo de demo

dataset = pd.read_csv(sys.argv[1])
path_output = sys.argv[2]

dataset = dataset[:36]

start_antigen = 0
start_antibody_L = 12
start_antibody_H = 24

matrix_fft_product = []

for i in range(len(dataset)):
    row_encoder = [dataset[value][i] for value in dataset.columns[:36]]
    antigen_row = np.array(row_encoder[start_antigen:start_antibody_L]).reshape(3,4)
    antibody_L = np.array(row_encoder[start_antibody_L: start_antibody_H]).reshape(4,3)
    antibody_H = np.array(row_encoder[start_antibody_H:]).reshape(3,4)

    matrix_antigen_antibodyL = antigen_row@antibody_L
    matrix_result = matrix_antigen_antibodyL@antibody_H

    vector_result = matrix_result.reshape(1, 12)
    matrix_fft_product.append(vector_result[0])

header = ["P_{}".format(i+1) for i in range(len(matrix_fft_product[0]))]
df_export = pd.DataFrame(matrix_fft_product, columns=header)
df_export.to_csv("{}export_dot_matrix_fft.csv".format(path_output), index=False)
    
