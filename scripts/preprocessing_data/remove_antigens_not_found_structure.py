import pandas as pd
import sys

print("Reading dataset")
dataset_actual = pd.read_csv(sys.argv[1]) #actual dataset with filter antigens by length
dataset_not_found = pd.read_csv(sys.argv[2])#antigens without structure
name_output = sys.argv[3]#name export dataset

dataset_actual = dataset_actual.loc[dataset_actual['filter_length'] == 1]
dataset_actual = dataset_actual.reset_index()

matrix_filter = []

print("Processing search")
for i in range(len(dataset_actual)):

    cont=0
    for j in range(len(dataset_not_found['id'])):
        if dataset_actual['id_seq'][i] == dataset_not_found['id'][j]:
            cont=1
            break 
    
    if cont==0:
        row = [dataset_actual['id_seq'][i], dataset_actual['seq'][i], 1]
        matrix_filter.append(row)

print("Export data")
df_export = pd.DataFrame(matrix_filter, columns=['id_seq', 'seq', 'filter_length'])
print(len(df_export))
df_export.to_csv(name_output, index=False)
