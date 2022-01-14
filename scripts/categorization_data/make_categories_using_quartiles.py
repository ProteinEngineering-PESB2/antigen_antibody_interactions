import pandas as pd 
import sys 
import numpy as np 

def get_class (value, ranges_data):
    if value < ranges_data['q1']:
        return 'C1'
    elif value >= ranges_data['q1'] and value<ranges_data['q2']:
        return 'C2'
    elif value >= ranges_data['q2'] and value<ranges_data['q3']:
        return 'C3'
    else:
        return 'C4'

def get_intervals (df):

    q1 = np.quantile(df['value_intensity'], .25)
    q2 = np.quantile(df['value_intensity'], .50)
    q3 = np.quantile(df['value_intensity'], .75)

    return [q1, q2, q3]

dataset = pd.read_csv(sys.argv[1])
data_export = sys.argv[2]

antibody_list = dataset['antibody'].unique()
range_data = {}

for antibody in antibody_list:

    df = dataset.loc[dataset['antibody'] == antibody]
    df = df.reset_index()
    row = get_intervals(df)
    range_data.update({antibody: {
                        "q1" : row[0],
                        "q2" : row[1],
                        "q3" : row[2],
                        } 
                    }) 

data_categories = []
for i in range(len(dataset)):
    category = get_class(dataset['value_intensity'][i], range_data[dataset['antibody'][i]])
    data_categories.append(category)

dataset['categories'] = data_categories
dataset.to_csv(data_export, index=False)
