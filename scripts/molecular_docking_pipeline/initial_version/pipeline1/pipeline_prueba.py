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

#based on blast results, select chain from antigen and save it into pdb file

def get_chain_antigen(path_structure_antigen, name_file, chain_antigen):
	cmd.load(path_structure_antigen+name_file)
	pdb_code=name_file.replace(".ent", "").replace("pdb","")
	cmd.select("chainAntigen", "pdb"+str(pdb_code)+" and chain "+chain_antigen)
	cmd.alter(("chain "+chain_antigen), "chain = 'A'")
	cmd.indicate("hetatm")
	cmd.remove("indicate")
	cmd.save(path_structure_antigen+"/"+pdb_code+"_A.pdb", "chainAntigen", -1, "pdb")

#makes folder for temporal files corresponding to complex... copies the antigen pdb file made before, the antibody heavy and light chain into the folder
#returns a list of information, which is used in the following functions

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

#for all the files in the folder made for the complex, it's implemented pdbfixer ... It removes non heteroatoms from structures (in case they aren't erased in the chain selection function)
#It also add the hidrogens to the structures, and complete possible missing atoms in residues.. It writes the corresponding repaires structures into pdb files

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

#The name of the chains for antibody light and heavy chain are changed with pymol ...  It also uses a base structure with both chains to attempt an initial conformation
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
			
#The antigen chain, antibody heavy and light chains are loaded and saved as one file in pymol ... this exports the first approach to the complex
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

#The complex generated before is renumbered, to abvoid problems in the repack and docking process, it renumbers from the 1 respect to the start in the pdb file
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

#outputs the command for repacking into the terminal
def repack_complex(path_results, structure_complex, repack_options, repack_protocol):
	
	os.system("srun ../../../../home/lmod/software/MPI/intel/2019.2.187-GCC-8.2.0-2.31.1/impi/2019.2.187/Rosetta/2019.07.60616/bin/rosetta_scripts.mpi.linuxiccrelease @"+repack_options+" -s "+structure_complex+" -parser:protocol "+repack_protocol+" -corrections::restore_talaris_behavior True -out:path:all "+path_results+" -nstruct 50")

#Takes the resulting scroting file and extract the best result based in the total score and the intra repulsion score. The lower the score, the best (is the sum of the energies)
#the intra repulsion refers to the energy of lennard jonnes between atoms in the same aminoacid
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

#Based on the cdr predictions for the antibodies, the aminoacids predicted are included in the protocol to avoid repacking of this aminoacids
#As the rest of the protocol don change, for any complex, just this section will be updated in base of the antibodies
def make_protocol_docking(cdr_prediction, path_results, structure_selected, docking_protocol):
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

#change the structure parameter in the file of docking options
def make_docking_options(docking_options, structure_selected):
	docking_file=open(docking_options, "r")
	list_of_lines = docking_file.readlines()
	list_of_lines[5] = "-s "+structure_selected+"\n"
	a_file = open(docking_options, "w")
	a_file.writelines(list_of_lines)
	a_file.close()

#Output the command for docking in the terminal
def run_docking_complex(path_results, structure_selected, docking_options, docking_protocol):
	os.system("srun ../../../../home/lmod/software/MPI/intel/2019.2.187-GCC-8.2.0-2.31.1/impi/2019.2.187/Rosetta/2019.07.60616/bin/rosetta_scripts.mpi.linuxiccrelease @"+docking_options+" -database ../../../../home/lmod//software/MPI/intel/2019.2.187-GCC-8.2.0-2.31.1/impi/2019.2.187/Rosetta/2019.07.60616/database/ -parser:protocol "+docking_protocol+" -out:suffix _full -corrections::restore_talaris_behavior True-out:path:all "+path_results+" -nstruct 1000")

#reads the scoring file generated for the docking, and selects the best complex based on the total score and the interface score ...
#Acording to rosetta, the best interface score are between -10 and -5
def select_best_docking(scoring_file_docking):
	docking_selected=""
	file_results=open(scoring_file, "r").read().splitlines()[2:]
	score=0
	i_sc=-5
	for line in file_results:
		results= re.sub(' +', ' ',line).split(" ")
		if float(results[1]) < score and float(results[5]) < rep_force and float(results[5]) > -10:
			docking_selected=results[-1]

	return docking_selected

#Move the selected docking to the directory to save the best results, the other docking structures, the folder of the complex and the repack structures are erased

def move_file_selected(path_results, docking_selected, path_dockings):
	shutil.copyfile(path_results+docking_selected, path_dockings+docking_selected)
	try:
		os.rmdir(path_complexes)
	except OSError as e:
		print("Error: %s : %s" % (dir_path, e.strerror))
	files = glob.glob(path_results+"*")
	for f in files:
		os.remove(f)


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
repack_options=parameters_file[9]
repack_protocol=parameters_file[10]
path_results=parameters_file[11]
cdr_prediction=pd.read_csv(parameters_file[12])
docking_options=parameters_file[13]
docking_protocol=parameters_file[14]
path_dockings=parameters_file[15].replace("\n","")


get_chain_antigen(path_structure_antigen, name_file, chain_antigen)
row_info=make_folder_and_copy_structures(path_structure_antibody, path_structure_antigen, interaction_antibody, interaction_antigen, info_antigens, path_complexes)
repair_pdbs_to_complex(row_info, path_complexes)
change_chains_antibodies(row_info, path_complexes, structure_align)
make_complex_pdb(row_info,path_complexes)
renumber_chains(path_complexes, row_info)
