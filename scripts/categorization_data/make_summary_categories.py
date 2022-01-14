import pandas as pd
import sys 
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots

dataset = pd.read_csv(sys.argv[1])
path_export = sys.argv[2]

unique_categories = dataset['categories'].unique()
counts = []

print("Make summary")

for i in range(len(unique_categories)):
    df = dataset.loc[dataset['categories'] == unique_categories[i]]
    counts.append(len(df))

fig = px.pie(dataset, values=counts, names=unique_categories)
fig.write_image(path_export+"summary_categories_using_q.png", scale=4, width=600, height=500)

print("Check each antibody")
antibody_list = dataset['antibody'].unique()

matrix_data = []

for antibody in antibody_list:
    df = dataset.loc[dataset['antibody'] == antibody]
    df = df.reset_index()

    unique_categories = df['categories'].unique()

    print("Make summary", antibody)

    for i in range(len(unique_categories)):
        df2 = df.loc[df['categories'] == unique_categories[i]]
        row = [antibody, unique_categories[i], len(df2)]
        matrix_data.append(row)
df = pd.DataFrame(matrix_data, columns=['antibody', 'category', 'counts'])

print(df)

fig2 = px.histogram(df, x="antibody", y="counts",
             color='category', barmode='group',
             height=400)
fig2.write_image(path_export+"summary_categories_per_antibody_using_q.png", scale=4, width=1600, height=800)

