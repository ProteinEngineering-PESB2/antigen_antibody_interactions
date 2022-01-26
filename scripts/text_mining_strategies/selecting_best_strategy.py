import os
import sys
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

def selecting_best_models(df, value, key):
    df_search = df.loc[df[key] > value]
    df_search.sort_values(by=['calinski', 'siluetas'])
    return df_search

print("Get files to process")
path_input = sys.argv[1]
dir_files = os.listdir(path_input)
list_files = [element for element in dir_files if "exploratory_clustering" in element]

list_df = []

for element in list_files:
    data_read = pd.read_csv(path_input+element)
    data_read['file'] = element
    list_df.append(data_read)

print("Join data")
data_join  = pd.concat(list_df)
data_join.reset_index(inplace=True)

matrix_data = []
for i in range(len(data_join)):
    if data_join['calinski'][i] != 'ERROR' and data_join['siluetas'][i] != 'ERROR':
        row = [data_join['algorithm'][i], data_join['params'][i], float(data_join['calinski'][i]), float(data_join['siluetas'][i]), data_join['file'][i]]
        matrix_data.append(row)

df_values = pd.DataFrame(matrix_data, columns=['algorithm', 'params', 'calinski', 'siluetas', 'file'])

print("Making histogram")
fig = make_subplots(rows=1, cols=2)
#['', 'percent', 'probability', 'density', 'probability density']
trace0 = go.Histogram(x=df_values['calinski'], name='Calinski-Harabasz', histnorm='probability')
trace1 = go.Histogram(x=df_values['siluetas'], name='Silhouettes coefficient', histnorm='probability')

fig.append_trace(trace0, 1, 1)
fig.append_trace(trace1, 1, 2)
fig.write_image(path_input+"histograms_results.svg")

print("Get max values")
min_cal = np.mean(df_values['calinski']) + 3*np.std(df_values['calinski'])
min_sil = np.mean(df_values['siluetas']) + 3*np.std(df_values['siluetas'])

best_cal = selecting_best_models(df_values, min_cal, 'calinski')
best_sil = selecting_best_models(df_values, min_sil, 'siluetas')

print("Export selected models")
best_sil.to_csv(path_input+"exporting_selected_best_silhouette.csv", index=False)
best_cal.to_csv(path_input+"exporting_selected_best_calinski.csv", index=False)