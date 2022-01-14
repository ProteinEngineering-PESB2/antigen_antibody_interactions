import pandas as pd
import sys
import plotly.graph_objects as go
import numpy as np 
import math 

dataset = pd.read_csv(sys.argv[1])
name_export = sys.argv[2]

list_distributions = {}

key_value = 'scaler_value'  #'value_intensity'
print("Get intensities")
antibody_list = dataset['antibody'].unique()
for antibody in antibody_list:

    df =  dataset.loc[dataset['antibody'] == antibody]
    values = [float(value) for value in df[key_value] if str(value) != "nan" and value is not None]
    list_distributions.update({antibody:values})

fig = go.Figure()

for antibody in antibody_list:
    fig.add_trace(go.Box(y=list_distributions[antibody], name=antibody))

fig.write_image(name_export, width=1600, height=950)
