#!/bin/bash

python3 pipeline_make_complexes A097-IOH28269.txt

srun ../../../../home/lmod/software/MPI/intel/2019.2.187-GCC-8.2.0-2.31.1/impi/2019.2.187/Rosetta/2019.07.60616/bin/rosetta_scripts.mpi.linuxiccrelease @repack.options -s structure_complex -parser:protocol repack.xml -corrections::restore_talaris_behavior True -out:path:all path_results -nstruct 50

python3 process_repack_and_make_docking path_results repack.fasc cdr3_antibodies.csv docking_full3.xml docking.options

srun ../../../../home/lmod/software/MPI/intel/2019.2.187-GCC-8.2.0-2.31.1/impi/2019.2.187/Rosetta/2019.07.60616/bin/rosetta_scripts.mpi.linuxiccrelease @docking.options -database ../../../../home/lmod//software/MPI/intel/2019.2.187-GCC-8.2.0-2.31.1/impi/2019.2.187/Rosetta/2019.07.60616/database/ -parser:protocol docking_full3.xml -out:suffix _full -corrections::restore_talaris_behavior True-out:path:all path_results -nstruct 10000

python3 select_docking_and_save dockinf_full.fasc path_results path_dockings

