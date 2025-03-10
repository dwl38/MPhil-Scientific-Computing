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

#python run.py --name='finetune_multihead' --model='../../finetune_multihead/model_finetune_multihead.model'
python run.py --name='mp0' --model='../../../WA-1/00-common/mace-models/2023-12-03-mace-128-L1_epoch-199.model'
python analyze.py mp0 --index='-1000:'

