import pandas as pd
import sys
import numpy as np

dataset = pd.read_csv(sys.argv[1])
dataset_antigens = pd.read_csv(sys.argv[2])
path_export = sys.argv[3]

print("Get filter sequences")
antigen_filters = dataset_antigens.loc[dataset_antigens['filter_length'] == 1]
antigen_filters.reset_index(inplace=True)

print("Get interactions using filter values")
antigens_filter_list = [id_value for id_value in antigen_filters['id_seq']]

filter_interaction = []
cont=0

for i in range(len(dataset)):
    if dataset['antigen'][i] in antigens_filter_list:
        filter_interaction.append(1)
        cont+=1
    else:
        filter_interaction.append(0)

dataset['filter'] = filter_interaction
dataset_f = dataset.loc[dataset['filter'] == 1]

antibody_list = dataset['antibody'].unique()

dict_summary = []

print("Searching data")
for antibody in antibody_list:

    df =  dataset_f.loc[dataset_f['antibody'] == antibody]
    min_value = np.min(df['value_intensity'])
    max_value = np.max(df['value_intensity'])
    mean_value = np.mean(df['value_intensity'])
    std_value = np.std(df['value_intensity'])
    q1 = np.quantile(df['value_intensity'], .25)
    q3 = np.quantile(df['value_intensity'], .75)

    row = [antibody, min_value, max_value, mean_value, std_value, q1, q3]
    dict_summary.append(row)

print("Export summary")
df_export = pd.DataFrame(dict_summary, columns=['antibody', 'min', 'max', 'mean', 'std', 'q1', 'q3'])
df_export.to_csv(path_export+"statistic_summary.csv", index=False)    

print("Export filter data")
dataset_f.to_csv(path_export+"filter_interactions.csv", index=False)