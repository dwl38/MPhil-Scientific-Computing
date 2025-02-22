#!/bin/bash

python -m venv ./
source ./bin/activate
sleep 5

pip install torch
pip install cuequivariance==0.1.0
pip install cuequivariance-torch==0.1.0
pip install cuequivariance-ops-torch-cu12==0.1.0

git clone https://github.com/ACEsuit/mace.git
pip install ./mace

