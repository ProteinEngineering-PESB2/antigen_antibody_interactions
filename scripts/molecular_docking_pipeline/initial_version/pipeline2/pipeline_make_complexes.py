import sys
import os
import pandas as pd
from pymol import cmd
import shutil
from pdbfixer import PDBFixer
from simtk.openmm.app import PDBFile
import re
from Bio import SeqIO
import xml.etree.ElementTree as ET

def get_chain_antigen(path_structure_antigen, name_file, chain_antigen):
	cmd.load(path_structure_antigen+name_file)
	pdb_code=name_file.replace(".ent", "").replace("pdb","")
	cmd.select("chainAntigen", "pdb"+str(pdb_code)+" and chain "+chain_antigen)
	cmd.alter(("chain "+chain_antigen), "chain = 'A'")
	cmd.indicate("hetatm")
	cmd.remove("indicate")
	cmd.save(path_structure_antigen+"/"+pdb_code+"_A.pdb", "chainAntigen", -1, "pdb")

def make_folder_and_copy_structures(path_structure_antibody,path_structure_antigen,interaction_antibody, interaction_antigen, info_antigens, path_complexes):
	files_antibodies=[ file for file in os.listdir(path_structure_antibody) if interaction_antibody in file]
	row=[]
	if len(files_antibodies) ==2:
		antigen_data= info_antigens.loc[info_antigens["Id_antigen"] == interaction_antigen].values.tolist()
		print(antigen_data)
		if len(antigen_data) == 1:
			name_dir_complex= interaction_antibody+"-"+antigen_data[0][0].split(":")[-1]
			print(name_dir_complex)
			name_file=antigen_data[0][2].replace("pdb", "").replace(".ent", "")+"_A.pdb"
			os.mkdir(path_complexes+name_dir_complex)
			shutil.copyfile(path_structure_antibody+files_antibodies[0], path_complexes+name_dir_complex+"/"+files_antibodies[0])
			shutil.copyfile(path_structure_antibody+files_antibodies[1], path_complexes+name_dir_complex+"/"+files_antibodies[1])
			shutil.copyfile(path_structure_antigen+"/"+name_file, path_complexes+name_dir_complex+"/"+name_file)
			row=[name_dir_complex, name_file, antigen_data[0][3], "A", files_antibodies[0], "A", "H", files_antibodies[1], "A", "L" ]
	return row

def repair_pdbs_to_complex(row_info, path_complexes):
	path_complex=path_complexes+row_info[0]+"/"
	data_antigen=row_info[1:4]
	data_antibody_H= row_info[4:7]
	data_antibody_L= row_info[7:11]
	row=[data_antigen, data_antibody_H, data_antibody_L]
	for lista in row:
		print(lista)
		fixer = PDBFixer(filename=path_complex+lista[0])
		name_file=lista[0].replace(".pdb","")
		fixer.findMissingResidues()
		fixer.findNonstandardResidues()
		fixer.replaceNonstandardResidues()
		fixer.removeHeterogens(True)
		fixer.findMissingAtoms()
		fixer.addMissingAtoms()
		fixer.addMissingHydrogens(7.0)
		PDBFile.writeFile(fixer.topology, fixer.positions, open(path_complex+name_file+"_repaired.pdb", 'w'))

def change_chains_antibodies(row_info, path_complexes, structure_align):
	path_complex=path_complexes+row_info[0]+"/"

	data_antibody_H= row_info[4:7]
	data_antibody_L= row_info[7:11]
	row_info=[data_antibody_H, data_antibody_L]
	#print(row_info)
	for lista in row_info:
		name_file_load=lista[0].replace(".pdb","_repaired.pdb")
		cmd.load(path_complex+name_file_load)
		name_file=name_file_load.replace(".pdb","")
		cmd.alter("(chain "+lista[1]+")","chain ='"+lista[2]+"'")
		cmd.load(structure_align)
		file_pattern=structure_align.split("/")[-1].replace(".pdb", "")
		cmd.align(name_file, file_pattern)
		cmd.drag(name_file)
		cmd.save(path_complex+lista[0].replace(".pdb","_repaired.pdb"), name_file, -1, "pdb")
		cmd.reinitialize()
			
def make_complex_pdb(row_info, path_complexes):
	name_file_antigen=row_info[1].replace(".pdb", "_repaired.pdb")
	name_file_H=row_info[4].replace(".pdb", "_repaired.pdb")
	name_file_L=row_info[7].replace(".pdb", "_repaired.pdb")

	cmd.load(path_complexes+row_info[0]+"/"+name_file_antigen)
	cmd.load(path_complexes+row_info[0]+"/"+name_file_H)
	cmd.load(path_complexes+row_info[0]+"/"+name_file_L)
	cmd.select("complex", "all")
	cmd.save(path_complexes+row_info[0]+"/"+row_info[0]+".pdb", "complex", -1, "pdb")
	cmd.reinitialize()

def renumber_chains(path_complexes, row_info):
	first_line=open(path_complexes+row_info[0]+"/"+row_info[0]+".pdb").readline().split(" ")
	list_numbers=[]
	list_nothing=[]
	for cosa in first_line:
		try:
			number=int(cosa)
			list_numbers.append(number)
		except:
			list_nothing.append(cosa)
	cmd.load(path_complexes+row_info[0]+"/"+row_info[0]+".pdb")

	cmd.select("chainA", "chain A")
	cmd.select("chainH", "chain H")
	cmd.select("chainL", "chain L")
	large_chain=len(cmd.get_model("chainA").get_residues())
	heavy_chain=len(cmd.get_model("chainH").get_residues())
	
	cmd.alter("chainA", "resi=str(int(resi)-"+str(int(list_numbers[-1])-1)+")")
	cmd.alter("chainH", "resi=str(int(resi)+"+str(large_chain)+")")
	cmd.alter("chainL", "resi=str(int(resi)+"+str(large_chain+heavy_chain)+")")
	cmd.sort()
	cmd.select("complex", "all")
	cmd.save(path_complexes+row_info[0]+"/"+row_info[0]+"_renumbered.pdb","complex")
	cmd.reinitialize()


parameters_file= open(sys.argv[1], "r").readline().split(" ")


path_structure_antigen=parameters_file[0]
name_file=parameters_file[1]
chain_antigen=parameters_file[2]
path_structure_antibody= parameters_file[3]
interaction_antigen=parameters_file[4]
interaction_antibody=parameters_file[5]
info_antigens=pd.read_csv(parameters_file[6])
path_complexes=parameters_file[7]
structure_align=parameters_file[8]

get_chain_antigen(path_structure_antigen, name_file, chain_antigen)
row_info=make_folder_and_copy_structures(path_structure_antibody, path_structure_antigen, interaction_antibody, interaction_antigen, info_antigens, path_complexes)
repair_pdbs_to_complex(row_info, path_complexes)
change_chains_antibodies(row_info, path_complexes, structure_align)
make_complex_pdb(row_info,path_complexes)
renumber_chains(path_complexes, row_info)
