import pandas as pd 
import sys 
import plotly.graph_objects as go

dataset = pd.read_csv(sys.argv[1])
export_file = sys.argv[2]

dataset = dataset.drop(columns=['antibody'])

matrix_data = []

for i in range(len(dataset)):
    row = [dataset[key][i] for key in dataset.columns]
    matrix_data.append(row)

fig = go.Figure(data=go.Heatmap(
                   z=matrix_data,
                   x=dataset.columns,
                   y=dataset.columns,
                   hoverongaps = False,
                   colorscale='Viridis'))
fig.write_image(export_file, width=1600, height=950)
