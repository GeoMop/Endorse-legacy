#!/bin/bash

# !! Needs redis server installed.
# sudo apt install redis

# download large data
tests/test_data/download.sh

# Setup virtual environment for development


python3 -m venv venv
source venv/bin/activate

python -m pip install -e submodules/bgem
python -m pip install -e submodules/bih
python -m pip install -e submodules/redis-simple-cache-3k
python -m pip install -r requirements.txt
