import pandas as pd
import sys

def get_class (value, ranges_data):
    if value >= ranges_data['min'] and value<ranges_data['q1']:
        return 'C1'
    elif value >= ranges_data['q1'] and value<ranges_data['mean']:
        return 'C2'
    elif value >= ranges_data['mean'] and value<ranges_data['q3']:
        return 'C3'
    else:
        return 'C4'

dataset = pd.read_csv(sys.argv[1])
summary_dataset = pd.read_csv(sys.argv[2])
data_export = sys.argv[3]

data_categories = []

print("Process ranges")
ranges = {}
for i in range(len(summary_dataset)):
    ranges.update({
                    summary_dataset['antibody'][i] :{
                        "min": summary_dataset['min'][i], 
                        "max":summary_dataset['max'][i],
                        "mean":summary_dataset['mean'][i],
                        "q1":summary_dataset['q1'][i],
                        "q3":summary_dataset['q3'][i]
                        }
                })
    
for i in range(len(dataset)):
    category = get_class(dataset['value_intensity'][i], ranges[dataset['antibody'][i]])
    data_categories.append(category)

dataset['categories'] = data_categories
dataset.to_csv(data_export, index=False)