#!/bin/bash

#SBATCH -J finetune-MACE
#SBATCH -A NIKIFORAKIS-LIM-SL2-GPU
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --time=36:00:00
#SBATCH --mail-type=NONE
##SBATCH --no-requeue
#SBATCH -p ampere

. /etc/profile.d/modules.sh                # Leave this line (enables the module command)
module purge                               # Removes all modules still loaded
module load rhel8/default-amp              # REQUIRED - loads the basic environment
module load cudnn

source ../venv/bin/activate

python ../venv/mace/mace/cli/run_train.py \
	--name="naive_finetuned" \
	--foundation_model="./mace_agnesi_medium.model" \
	--multiheads_finetuning=False \
	--train_file="../data/gra-wat-protonic-defects.xyz" \
	--valid_fraction=0.05 \
	--energy_weight=1.0 \
	--forces_weight=1.0 \
	--E0s="{1:-13.656315639066156,6:-146.65071522867848,8:-432.2257363227328}" \
	--model="MACE" \
	--lr=0.01 \
	--scaling="rms_forces_scaling" \
	--batch_size=10 \
	--max_num_epochs=1500 \
	--ema \
	--ema_decay=0.99 \
	--amsgrad \
	--device=cuda


