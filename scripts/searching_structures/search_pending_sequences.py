import pandas as pd 
import os 
import sys 

data_with_all_seqs = pd.read_csv(sys.argv[1])
path_to_process_sequences = sys.argv[2]
path_export = sys.argv[3]

list_files = os.listdir(path_to_process_sequences)

list_sequences_processed = []

print("Read processed sequences")
for element in list_files:
    df = pd.read_csv(path_to_process_sequences+element)
    id_seq = df['Id_sequence'][0]
    list_sequences_processed.append(id_seq)

print("searching pending sequences")
matrix_data = []

for i in range(len(data_with_all_seqs)):
    if data_with_all_seqs['id_seq'][i] not in list_sequences_processed and data_with_all_seqs['filter_length'][i] == 1:
        row = [data_with_all_seqs['id_seq'][i], data_with_all_seqs['seq'][i]]
        matrix_data.append(row)

df_pending = pd.DataFrame(matrix_data, columns=['id_seq', 'seq'])

print("Export pending sequnces")
df_pending.to_csv(path_export+"pending_sequences.csv", index=False)