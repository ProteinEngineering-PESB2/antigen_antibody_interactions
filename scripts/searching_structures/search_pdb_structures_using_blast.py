import sys
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import pandas as pd
import multiprocessing as mp
import numpy as np

def search_blast_sequence(id_sequence, sequence, name_output, header):
    print("Process sequence: ", id_sequence)
    matrix_export = []
    try:
        result_handle = NCBIWWW.qblast("blastp", "pdb", sequence)
        blast_records = NCBIXML.read(result_handle)
        
        for alignment in blast_records.alignments:
            for hsp in alignment.hsps:
                if "Homo sapiens" in alignment.title and hsp.expect <= 0.5:	
                    result= alignment.title.split("pdb")
                    for cosa in result:
                        if len(cosa) > 1:
                            row = []
                            row.append(id_sequence)
                            row.append(cosa.split("|")[1])
                            row.append(cosa.split("|")[2].split(" ")[0])
                            row.append(alignment.length)
                            row.append(hsp.expect)
                            matrix_export.append(row)
    except:
        pass

    df_export = pd.DataFrame(matrix_export, columns=header)
    df_export.to_csv(name_output, index=False)

def process_element_in_range(antigen_sequences, start, stop, path_output):
    header=["Id_sequence", "Id_PDB", "Chain", "Alignment_length", "evalue"]

    for i in range(start, stop):
        name_output = "{}{}.csv".format(path_output, antigen_sequences['id_seq'][i])
        name_output = name_output.replace(":", "_").replace("~", "_")
        sequence = antigen_sequences['seq'][i]
        id_sequence = antigen_sequences['id_seq'][i]

        search_blast_sequence(id_sequence, sequence, name_output, header)

antigen_sequences= pd.read_csv(sys.argv[1])
path_output=sys.argv[2]

cpu_number= 50
#we process the script using paralel option
number_data = int(len(antigen_sequences)/cpu_number)
rest_element = int(len(antigen_sequences)%cpu_number)

#get index
index_data = []

print("Evaluate split process")
for i in range(cpu_number):
    start = i*number_data
    stop = start+number_data

    if i == cpu_number-1:
        stop = stop+rest_element
    row = [start, stop]
    index_data.append(row)

#star paralelism
with mp.Manager() as manager:
    
    print("Init process")
    
    process_sequences = []
    for i in range(cpu_number):

        #create vector list (for 1)
        start = index_data[i][0]
        stop = index_data[i][1]
        print(start, stop)

        #paralelize process
        process_sequences.append(mp.Process(target=process_element_in_range, args=[antigen_sequences, start, stop, path_output]))

    #run and join data process	
    for p in process_sequences:
        print("Start process evaluate sequences")
        p.start()

    for p in process_sequences:
        print("Join data evaluate sequences process")
        p.join()
        