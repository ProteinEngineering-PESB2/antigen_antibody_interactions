import pandas as pd
import sys
import numpy as np 
import os
import json

def get_statistic_values(dataset):

    dict_values = {}

    #get points
    points = [value for value in dataset.columns if "P_" in value]
    dict_values.update({"points": points})

    #get statistics
    for point in points:
        mean_value = np.mean(dataset[point])
        std_value = np.std(dataset[point])
        min_value = np.min(dataset[point])
        max_value = np.max(dataset[point])
        low_IC = mean_value - 1.96*std_value
        high_IC = mean_value + 1.96*std_value
        iq_low = np.quantile(dataset[point], .25)
        iq_high = np.quantile(dataset[point], .75)

        point_dict = {"mean": mean_value, "std_data": std_value, "min_value": min_value, "max_value": max_value, "low_ic" : low_IC, "high_ic" : high_IC, "low_iq" : iq_low, "high_iq" : iq_high}
        dict_values.update({point:point_dict})

    return dict_values

path_input = sys.argv[1]
path_output = sys.argv[2]

list_datasets = os.listdir(path_input)

dict_response = {}

for file_process in list_datasets:
    name_dataset = "{}{}".format(path_input, file_process)
    print("Process ", name_dataset)
    dataset = pd.read_csv(name_dataset)
    
    response_statistics = get_statistic_values(dataset)
    dict_response.update({file_process.split(".")[0]: response_statistics})

print("Export json")
with open('{}summary_json_values.json'.format(path_output), 'w') as outfile:
    json.dump(dict_response, outfile)




