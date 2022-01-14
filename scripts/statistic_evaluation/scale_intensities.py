import pandas as pd
import sys 
from sklearn.preprocessing import MinMaxScaler
import numpy as np 

def scale_data (data):
    min_value = np.min(data)
    max_value = np.max(data)
    scaler_values = []

    for x in data:
        X_std = (x - min_value) / (max_value - min_value)    
        scaler_values.append(X_std)
    
    return scaler_values

def scale_data_using_zscore(data):
    mean = np.mean(data)
    std = np.std(data)
    scaler_values = []

    for x in data:
        X_std = (x - mean) / std
        scaler_values.append(X_std)
    
    return scaler_values

dataset = pd.read_csv(sys.argv[1])
path_export = sys.argv[2]

list_distributions = {}

print("Get intensities")
antibody_list = dataset['antibody'].unique()
for antibody in antibody_list:

    df =  dataset.loc[dataset['antibody'] == antibody]
    df.reset_index(inplace=True)
    values = [float(value) for value in df['value_intensity'] if str(value) != "nan" and value is not None]
    interactions = [df['antigen'][i] for i in range(len(df)) if str(df['antigen'][i])!= "nan" and str(df['antigen'][i]) is not None]
    list_distributions.update({antibody:{"values":values, "antigens":interactions}})

print("Start scaler")
matrix_export = []
for antibody in antibody_list:
    print("Scaler antibody values", antibody)
    
    scaler_values = scale_data_using_zscore(list_distributions[antibody]['values'])
    antigens = list_distributions[antibody]['antigens']
    
    for i in range(len(scaler_values)):
        row = [antibody, antigens[i], scaler_values[i]]
        matrix_export.append(row)

df_export = pd.DataFrame(matrix_export, columns=['antibody', 'antigen', 'scaler_value'])
df_export.to_csv(path_export+"scaler_interactions_zscore.csv", index=False)