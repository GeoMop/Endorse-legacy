#!/bin/bash

# !! Needs redis server installed.
# sudo apt install redis

# Setup virtual environment for development

python -m venv venv
source venv/bin/activate

python -m pip install -e ../bgem
python -m pip install -e ../bih
python -m pip install -e ../redis-simple-cache-3k
python -m pip install -r requirements.txt
