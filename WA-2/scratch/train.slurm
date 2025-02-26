#!/bin/bash

#SBATCH -J scratch-MACE
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
	--name="model_from_scratch" \
	--train_file="../data/gra-wat-protonic-defects.xyz" \
	--valid_fraction=0.05 \
	--config_type_weights='{"Default":1.0}' \
	--E0s="{1:-13.656315639066156,6:-146.65071522867848,8:-432.2257363227328}" \
	--model="MACE" \
	--hidden_irreps='128x0e + 128x1o' \
	--r_max=6.0 \
	--batch_size=10 \
	--max_num_epochs=1500 \
	--swa \
	--start_swa=1200 \
	--ema \
	--ema_decay=0.99 \
	--amsgrad \
	--restart-latest \
	--device=cuda


