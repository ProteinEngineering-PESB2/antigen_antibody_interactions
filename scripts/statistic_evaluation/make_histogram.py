import pandas as pd
import sys
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np 
import math 

dataset = pd.read_csv(sys.argv[1])
name_export = sys.argv[2]

list_distributions = {}

key_value =  'value_intensity' #'scaler_value'  #'value_intensity'
print("Get intensities")
antibody_list = dataset['antibody'].unique()
for antibody in antibody_list:

    df =  dataset.loc[dataset['antibody'] == antibody]
    values = [float(value) for value in df[key_value] if str(value) != "nan" and value is not None]
    list_distributions.update({antibody:values})

fig = make_subplots(rows=5, cols=9)

list_traces = []

for antibody in antibody_list:
    trace = go.Histogram(x=list_distributions[antibody], autobinx=True)
    list_traces.append(trace)

index=0
for i in range(1, 6):
    for j in range(1, 10):
        fig.append_trace(list_traces[index], i, j)
        index+=1

fig.write_image(name_export, width=1600, height=950)
