#!/bin/bash

#SBATCH -J crystal-rtp
#SBATCH -A NIKIFORAKIS-LIM-SL2-GPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --time=24:00:00
#SBATCH --mail-type=NONE

#SBATCH -p ampere

. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel8/default-amp              # REQUIRED - loads the basic environment

source ../00-common/venv/bin/activate

python run.py

