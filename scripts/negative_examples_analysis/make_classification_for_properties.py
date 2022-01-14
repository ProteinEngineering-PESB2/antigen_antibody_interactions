import pandas as pd
import sys
import json
import os

def get_pbb_by_property(dataset, config_pbb, property_value):

    points = config_pbb[property_value]['points']

    matrix_data = []

    for i in range(len(dataset)):
        count = 0
        id_seq = dataset['id_sequence'][i]

        for point in points:

            value_point = dataset[point][i]
            ic_high = config_pbb[property_value][point]['high_ic']
            ic_low = config_pbb[property_value][point]['low_ic']

            if value_point <= ic_high and value_point>= ic_low:
                count+=1

        pbb_value = float(count)/float(len(points))
        row = [id_seq, pbb_value]
        matrix_data.append(row)
    
    df_data = pd.DataFrame(matrix_data, columns=['id_seq', 'pbb_{}'.format(property_value)])
    return df_data

path_encoding_to_evaluate = sys.argv[1]
config_json = sys.argv[2]
path_to_export = sys.argv[3]

#read summary data in json file
with open(config_json) as json_data:
    data_config = json.load(json_data)

probability_matrix = []

#csvs to process
list_element = os.listdir(path_encoding_to_evaluate)

for element in list_element:

    name_dataset = "{}{}".format(path_encoding_to_evaluate, element)
    print("Process ", name_dataset)
    dataset = pd.read_csv(name_dataset)
    
    property_value = element.split(".")[0]
    df_response = get_pbb_by_property(dataset, data_config, property_value)

    probability_matrix.append(df_response)

df_merge = probability_matrix[0]

for i in range(1, len(probability_matrix)):
    df_merge = pd.merge(df_merge, probability_matrix[i], on='id_seq')

#get full pbb
full_pbb = []

for i in range(len(df_merge)):

    id_seq = df_merge['id_seq'][i]
    pbb_value = 1

    for column in df_merge.columns:
        if column != 'id_seq':
            pbb_value = pbb_value*df_merge[column][i]
    
    row = [id_seq, pbb_value]
    full_pbb.append(row)

df_export = pd.DataFrame(full_pbb, columns=['id_seq', 'pbb_full'])
df_export.to_csv("{}estimated_latent_pbb.csv".format(path_to_export), index=False)
