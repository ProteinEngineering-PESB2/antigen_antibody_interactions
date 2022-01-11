import pandas as pd
from multiprocessing import Pool
import os
from Bio import SeqIO, Seq
from Bio.SeqRecord import SeqRecord
import sys
chromedriver_path = "/usr/bin/chromedriver"
import re

def search_blast_sequence(id_sequence, sequence):
    sequence = sequence.replace("*","")
    record = SeqRecord(Seq.Seq(sequence), id_sequence)
    file = "fasta/{}.fasta".format(id_sequence)
    SeqIO.write(record, file, "fasta")
    output = "uniprot/{}.result".format(id_sequence)
    command= "blastp -db swissprot -query {} -evalue 0.5 -taxids 9606 -out {} ".format(file, output)
    os.system(command)
    f = open(output, "r")
    blastp_text = f.read()

    f.close()
    found = False
    if("No hits found" not in blastp_text):
        inicio = blastp_text.find("Value") + 7
        fin = blastp_text.find("\n>")
        tabla = blastp_text[inicio:fin]

        final_csv = "Sequence_id,Uniprot_id,Score,e-value\n"
        for row in tabla.split("\n"):
            if(row != ''):
                rex = re.compile(r'\s+')
                res = rex.sub(' ', row)
                res = res.strip().split(" ")
                final_csv += "{},{},{},{}\n".format(id_sequence,res[0], res[-2], res[-1])
                found = True
        f = open(output, "w")
        f.write(final_csv)
        f.close()
    else:
        os.system("rm {}".format(output))
    return found
    
def process(row):
    row = row[1]
    found = search_blast_sequence(row.id, row.sequence)
    
if __name__ == '__main__':
    dataset = sys.argv[1]
    try:
        os.mkdir("fasta")
    except:
        pass
    try:
        os.mkdir("uniprot")
    except:
        pass 
    try:
        os.mkdir("alphafold_pdb")
    except:
        pass
    dataset = pd.read_csv(dataset)
    with Pool(20) as p:
        p.map(process, dataset.iterrows())