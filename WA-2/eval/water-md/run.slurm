#!/bin/bash

#SBATCH -J water-rdf
#SBATCH -A NIKIFORAKIS-LIM-SL2-GPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --time=12:00:00
#SBATCH --mail-type=NONE

#SBATCH -p ampere

. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel8/default-amp              # REQUIRED - loads the basic environment

source ../../venv/bin/activate

python run.py --name='finetune_multihead' --model='../../finetune_multihead/model_finetune_multihead.model'
python analyze.py finetune_multihead --index='-1000:'

