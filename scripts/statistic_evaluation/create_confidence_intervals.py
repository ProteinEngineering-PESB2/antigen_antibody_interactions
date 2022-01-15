import scipy.stats as st
import numpy as np 
import pandas as pd 
import sys 

dataset = pd.read_csv(sys.argv[1])
name_export = sys.argv[2]

antibody_list = dataset['antibody'].unique()
key_value = 'scaler_value'  #'value_intensity'

dict_summary = []

print("Searching data")
for antibody in antibody_list:

    df =  dataset.loc[dataset['antibody'] == antibody]
    interval = st.t.interval(alpha=0.95, df=len(df)-1, loc=np.mean(df[key_value]), scale=st.sem(df[key_value])) 
    row = [antibody, interval[0], interval[1]]
    dict_summary.append(row)

print("Export summary")
df_export = pd.DataFrame(dict_summary, columns=['antibody', 'IC_low', 'IC_High'])
df_export.to_csv(name_export, index=False)    

