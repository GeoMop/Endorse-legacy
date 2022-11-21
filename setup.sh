#!/bin/bash

# !! Needs redis server installed.
# sudo apt install redis

# download large data
tests/test_data/download.sh

# Setup virtual environment for development


python3 -m venv venv
# python3 -m venv --system-site-packages venv
source venv/bin/activate

python -m pip install -e ../bgem
python -m pip install -e ../bih
python -m pip install -e ../redis-simple-cache-3k
python -m pip install -r requirements.txt
