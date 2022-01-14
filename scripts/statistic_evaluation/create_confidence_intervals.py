import scipy.stats as st
import numpy as np 
import pandas as pd 
import sys 

dataset = pd.read_csv(sys.argv[1])
path_export = sys.argv[2]

antibody_list = dataset['antibody'].unique()

dict_summary = []

print("Searching data")
for antibody in antibody_list:

    df =  dataset.loc[dataset['antibody'] == antibody]
    interval = st.t.interval(alpha=0.95, df=len(df)-1, loc=np.mean(df['value_intensity']), scale=st.sem(df['value_intensity'])) 
    row = [antibody, interval[0], interval[1]]
    dict_summary.append(row)

print("Export summary")
df_export = pd.DataFrame(dict_summary, columns=['antibody', 'IC_low', 'IC_High'])
df_export.to_csv(path_export+"confidence_intervals_summary.csv", index=False)    

