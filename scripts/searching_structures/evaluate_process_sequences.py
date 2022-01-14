import os
import sys
import pandas as pd

list_files = os.listdir(sys.argv[1])
list_files = [value for value in list_files if ".csv" in value]

elements_with_structure = []

print("Search element")
for element in list_files:

    dataset = pd.read_csv(sys.argv[1]+element)
    if len(dataset)>0:
        elements_with_structure.append(element)
print(len(elements_with_structure))

print("Move element with not structure")
for element in elements_with_structure:
    command = "mv {}{} ../results/new_results".format(sys.argv[1], element)
    os.system(command)