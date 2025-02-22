#!/bin/bash

cd ..
mkdir ./tmp_clean
mv ./venv/setup.sh ./tmp_clean/setup.sh
mv ./venv/clean.sh ./tmp_clean/clean.sh
rm -rf ./venv
mkdir ./venv
mv ./tmp_clean/setup.sh ./venv/setup.sh
mv ./tmp_clean/clean.sh ./venv/clean.sh
rm -rf ./tmp_clean
cd ./venv

