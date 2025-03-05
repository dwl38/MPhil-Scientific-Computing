
python ../../venv/mace/mace/cli/eval_configs.py \
    --configs="../../data/nacl-water-RPA.xyz" \
    --model="../../scratch/model_scratch.model" \
    --output="./scratch.xyz"

python ../../venv/mace/mace/cli/eval_configs.py \
    --configs="../../data/nacl-water-RPA.xyz" \
    --model="../../finetune_naive/model_finetune_naive.model" \
    --output="./finetune_naive.xyz"

python ../../venv/mace/mace/cli/eval_configs.py \
    --configs="../../data/nacl-water-RPA.xyz" \
    --model="../../finetune_multihead/model_finetune_multihead.model" \
    --output="./finetune_multihead.xyz"
