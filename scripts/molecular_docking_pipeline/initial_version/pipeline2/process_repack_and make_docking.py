import sys
import os
import pandas as pd
import re
from Bio import SeqIO
import xml.etree.ElementTree as ET

def process_results_repack(path_results, scoring_file_repack):
	structure_selected=""
	file_results=open(scoring_file, "r").read().splitlines()[2:]
	score=0
	rep_force=20
	for line in file_results:
		results= re.sub(' +', ' ',line).split(" ")
		if float(results[1]) < score and float(results[6]) < rep_force:
			structure_selected=results[-1]

	return path_results+structure_selected

def make_protocol_docking(cdr_prediction, structure_selected, docking_protocol):
	antibody=structure_selected.split("/")[-1].split("-")[0]
	row_loc=cdr_prediction.loc[cdr_prediction["id"] == antibody]
	seq=""
	for record in SeqIO.parse(structure_selected, "pdb-atom"):
		seq=seq+record.seq
	print(row_loc["cdr3_heavy"].values[0])
	heavy_cdr= seq.find(row_loc["cdr3_heavy"].values[0])+1
	light_cdr=seq.find(row_loc["cdr3_light"].values[0])+1
	list_indexes=[]
	for n in range(len(list(row_loc["cdr3_heavy"].values[0]))):
		list_indexes.append(str(heavy_cdr+n))
	for n in range(len(list(row_loc["cdr3_light"].values[0]))):
		list_indexes.append(str(light_cdr+n))
	mytree = ET.parse(docking_protocol)
	myroot = mytree.getroot()
	indexes_protocol=",".join(list_indexes)
	for neighbor in myroot.iter('PreventResiduesFromRepacking'):
		neighbor.set("residues", indexes_protocol)
		print(neighbor.attrib)
	mytree.write(docking_protocol)

def make_docking_options(docking_options, structure_selected):
	docking_file=open(docking_options, "r")
	list_of_lines = docking_file.readlines()
	list_of_lines[5] = "-s "+structure_selected+"\n"
	a_file = open(docking_options, "w")
	a_file.writelines(list_of_lines)
	a_file.close()

path_results=sys.argv[1]
scoring_file_repack=sys.argv[2]
cdr_prediction=pd.read_csv(sys.argv[3])
docking_protocol=sys.argv[4]
docking_options=sys.argv[5]

structure_selected=process_results_repack(path_results, scoring_file_repack)
make_protocol_docking(cdr_prediction, structure_selected, docking_protocol)
make_docking_options(docking_options, structure_selected)

#write to file the selected structure and the parameters for 