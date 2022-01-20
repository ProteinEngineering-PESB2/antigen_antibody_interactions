#!/bin/bash
#SBATCH --job-name=demo_docking
#SBATCH --partition=slims
#SBATCH --output=demo_docking_yasna_%j.out
#SBATCH --error=demo_docking_yasna_%j.err
#SBATCH --mail-user=ybarrera11@alumnos.utalca.cl
#SBATCH --mail-type=ALL

./run_demo2.sh
