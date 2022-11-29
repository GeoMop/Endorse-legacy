#!/bin/bash

# !! Needs redis server installed.
# sudo apt install redis

# download large data
(cd tests/test_data/ ; sh download.sh)

# Setup virtual environment for development


python3 -m venv venv
# python3 -m venv --system-site-packages venv
source venv/bin/activate

python -m pip install -e submodules/bgem
python -m pip install -e submodules/bih
python -m pip install -e submodules/redis-cache
python -m pip install -e submodules/MLMC
python -m pip install -r requirements.txt
