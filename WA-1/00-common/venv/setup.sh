#!/bin/bash

python -m venv ./
source ./bin/activate
sleep 5

pip install torch==2.6.0
pip install cuequivariance==0.2.0
pip install cuequivariance-torch==0.2.0
pip install cuequivariance-ops-torch-cu12==0.2.0

git clone https://github.com/ACEsuit/mace.git
pip install ./mace

