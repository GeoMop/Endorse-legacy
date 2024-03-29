#!/bin/bash
echo "Make sure you called 'source load_modules.sh' before!"
echo "Creating python environment.."
python3 -m venv env

#source ./load_modules.sh
source env/bin/activate
python --version
which python
which pip
pip install --upgrade pip
pip -V

pip install pyyaml attrs numpy pandas scipy
pip install git+https://github.com/GeoMop/bgem.git@JB_endorse_mlmc

#pip freeze
deactivate

