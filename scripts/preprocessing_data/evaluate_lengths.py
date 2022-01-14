import pandas as pd
import sys
import numpy as np 

print("Read dataset")
dataset = pd.read_csv(sys.argv[1])
path_export = sys.argv[2]

print("Estimated length")
dataset['length'] = [len(element) for element in dataset['seq']]

print("Get statistic values")
min_value = np.min(dataset['length'])
max_value = np.max(dataset['length'])
mean_value = np.mean(dataset['length'])
std_value = np.std(dataset['length'])
q1_value = np.quantile(dataset['length'], .25)
q3_value = np.quantile(dataset['length'], .75)

print(min_value, max_value)
print(mean_value, std_value)
print(q1_value, q3_value)

print("Define statistic values")
#min_value_f = mean_value - 1.5*std_value
max_value_f = mean_value + 1.5*std_value

#iqr = q3_value - q1_value
#min_value_f = q1_value - 1.5*iqr
#max_value_f = q3_value + 1.5*iqr

min_value_f = 50#largo minimo de una proteina
#max_value_f = 600

print("Define range")
print(min_value_f, max_value_f)

cont=0

filter_length = []
for i in range(len(dataset)):
    if dataset['length'][i]<= max_value_f and dataset['length'][i] >=min_value_f:
        filter_length.append(1)
    else:
        filter_length.append(0)

dataset['filter_length'] = filter_length
dataset.to_csv(path_export+"filter_antigens_by_length.csv", index=False)


