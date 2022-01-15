from scipy.stats import ks_2samp
import pandas as pd
import numpy as np 
import math
import sys 

def apply_test(data1, data2, significance_level):
    U1, p = ks_2samp(data1, data2)
    if p <= significance_level:
        return 1
    else:
        return 0

dataset = pd.read_csv(sys.argv[1])
name_export = sys.argv[2]

list_distributions = {}

key_value =  'value_intensity' #'scaler_value'  #'value_intensity'

print("Get intensities")
antibody_list = dataset['antibody'].unique()
for antibody in antibody_list:

    df =  dataset.loc[dataset['antibody'] == antibody]
    values = [float(value) for value in df[key_value] if str(value) != "nan" and value is not None]
    list_distributions.update({antibody:values})

matrix_comparison = []

print("Estimated bonferroni corrections")
alpha = 0.01
number_experiments = math.factorial(len(antibody_list))/(math.factorial(len(antibody_list)-2)*2)
significance_level = alpha/number_experiments

print(significance_level)

print("Compare distributions")

for antibody1 in antibody_list:
    row = []
    for antibody2 in antibody_list:
        if antibody1 == antibody2:
            row.append(0)
        else:
            response = apply_test(list_distributions[antibody1], list_distributions[antibody2], significance_level)
            row.append(response)
    matrix_comparison.append(row)
     

df_export = pd.DataFrame(matrix_comparison, columns=antibody_list)
df_export['antibody'] = antibody_list
df_export.to_csv(name_export, index=False)