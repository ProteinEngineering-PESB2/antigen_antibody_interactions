import pandas as pd
import sys

def get_class (value, ranges_data):
    if value < ranges_data['min']:
        return 'C1'
    elif value >= ranges_data['min'] and value<=ranges_data['max']:
        return 'C2'
    else:
        return 'C3'

dataset = pd.read_csv(sys.argv[1])
confidence_intervals = pd.read_csv(sys.argv[2])
data_export = sys.argv[3]

data_categories = []

print("Process ranges")
ranges = {}
for i in range(len(confidence_intervals)):
    ranges.update({confidence_intervals['antibody'][i] : {"min": confidence_intervals['IC_low'][i], "max":confidence_intervals['IC_High'][i]}})
    
for i in range(len(dataset)):
    category = get_class(dataset['value_intensity'][i], ranges[dataset['antibody'][i]])
    data_categories.append(category)

dataset['categories'] = data_categories
dataset.to_csv(data_export, index=False)