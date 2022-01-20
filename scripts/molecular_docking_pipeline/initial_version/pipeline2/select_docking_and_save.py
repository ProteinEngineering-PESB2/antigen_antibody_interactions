import sys
import os
import shutil
import re

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


def move_file_selected(path_results, docking_selected, path_dockings):
	shutil.copyfile(path_results+docking_selected, path_dockings+docking_selected)
	try:
		os.rmdir(path_complexes)
	except OSError as e:
		print("Error: %s : %s" % (dir_path, e.strerror))
	files = glob.glob(path_results+"*")
	for f in files:
		os.remove(f)


scoring_file_docking=sys.argv[1]
path_results=sys.argv[2]
path_dockings=sys.argv[3]

docking_selected=select_best_docking(scoring_file_docking)
move_file_selected(path_results, docking_selected, path_dockings)
