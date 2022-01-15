import pandas as pd
import sys
import numpy as np

dataset = pd.read_csv(sys.argv[1])
name_export = sys.argv[2]

print("Get filter sequences")
antibody_list = dataset['antibody'].unique()

dict_summary = []

key_value = 'scaler_value'  #'value_intensity'

print("Searching data")
for antibody in antibody_list:

    df =  dataset.loc[dataset['antibody'] == antibody]
    min_value = np.min(df[key_value])
    max_value = np.max(df[key_value])
    mean_value = np.mean(df[key_value])
    std_value = np.std(df[key_value])
    q1 = np.quantile(df[key_value], .25)
    q3 = np.quantile(df[key_value], .75)

    row = [antibody, min_value, max_value, mean_value, std_value, q1, q3]
    dict_summary.append(row)

print("Export summary")
df_export = pd.DataFrame(dict_summary, columns=['antibody', 'min', 'max', 'mean', 'std', 'q1', 'q3'])
df_export.to_csv(name_export, index=False)    
