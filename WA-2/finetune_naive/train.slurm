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
	--name="model_finetune_naive" \
	--train_file="../data/nacl-water-RPA.xyz" \
	--valid_fraction=0.05 \
	--config_type_weights='{"Default":1.0}' \
	--model="MACE" \
	--foundation_model="medium" \
	--multiheads_finetuning=False \
	--batch_size=4 \
	--max_num_epochs=100 \
	--ema \
	--ema_decay=0.995 \
	--amsgrad \
	--device=cuda \
	--keep_checkpoints \
	--save_all_checkpoints


