import pandas as pd
import os
import numpy as np
import sys
dataset = sys.argv[1]
sequences = pd.read_csv(dataset)
#print(sequences)
uniprot_folder = "./uniprot/"
uniprot_files = [file.replace(".result", "") for file in os.listdir(uniprot_folder)]
#print(uniprot_files)
structures = os.listdir("alphafold_pdb")

for row in sequences.iterrows():
    id = row[0]
    row = row[1]
    if(row.id in uniprot_files):
        try:
            result = pd.read_csv("./uniprot/{}.result".format(row.id), sep=",")
            sequences.loc[ id,"Has_Uniprot"] = True
            uniprot_id = result.loc[0].Uniprot_id[:-2]
            e_value = result.loc[0]["e-value"]
            score = result.loc[0]["Score"]
            sequences.loc[ id,"Uniprot_id"] = uniprot_id
            sequences.loc[ id,"Score"] = score
            sequences.loc[ id,"e-value"] = e_value
            structures = os.listdir("./alphafold_pdb")
            if(uniprot_id+".pdb" in structures):
                sequences.loc[ id,"Has_Alphafold_Structure"] = True

                f = open("alphafold_pdb/" + row.Uniprot_id + ".pdb")
                text = f.read()
                description = text.split("COMPND")[2].strip()[12:-1]
                sequences.loc[id,"Description"]= description
            else:
                sequences.loc[ id,"Has_Alphafold_Structure"] = False
        except:
            sequences.loc[ id,"Has_Uniprot"] = False
            sequences.loc[ id,"Uniprot_id"] = np.nan
            sequences.loc[ id,"Has_Alphafold_Structure"] = False
    else:
        sequences.loc[ id,"Has_Uniprot"] = False
        sequences.loc[ id,"Uniprot_id"] = np.nan
        sequences.loc[ id,"Has_Alphafold_Structure"] = False
        

sequences.to_csv("Pending_Sequences_Structures.csv", index=False)
sequences.query("Has_Alphafold_Structure == False").to_csv("Not_Found_Sequences.csv", index=False)

