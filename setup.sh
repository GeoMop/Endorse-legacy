#!/bin/bash

# !! Needs redis server installed.
# sudo apt install redis

# download large data
(cd tests/test_data/ ; sh download.sh)

# Setup virtual environment for development


python3 -m venv venv
#python3 -m venv --system-site-packages venv
source venv/bin/activate
python3 -m pip install -r requirements.txt

python3 -m pip install -e submodules/bgem
python3 -m pip install -e submodules/bih
python3 -m pip install -e submodules/redis-cache
python3 -m pip install -e submodules/MLMC

python3 -m pip install -e .
