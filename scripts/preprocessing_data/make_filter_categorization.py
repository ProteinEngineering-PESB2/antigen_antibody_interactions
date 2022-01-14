import pandas as pd
import sys

#input data
print("Reading data")
dataset = pd.read_csv(sys.argv[1])
dataset_antigens = pd.read_csv(sys.argv[2])
name_export = sys.argv[3]

print("Get filter sequences")
antigen_filters = dataset_antigens.loc[dataset_antigens['filter_length'] == 1]
print(len(antigen_filters))
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

print("Export dataset")
dataset_f.to_csv(name_export, index=False)