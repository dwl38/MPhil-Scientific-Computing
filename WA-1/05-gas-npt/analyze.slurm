#!/bin/bash

#SBATCH -J mol-count
#SBATCH -A SCHRAN-SL3-CPU
#SBATCH -p icelake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:45:00
#SBATCH --mail-type=NONE

. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel8/default-icl              # REQUIRED - loads the basic environment

conda init
source ~/.bashrc
conda activate venv

python ../00-common/analyze.py 1100K --filter TATB TATB-F1 TATB-F2 TATB-F3 H2O NO2


