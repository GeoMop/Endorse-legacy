#!/bin/bash

# Setup virtual environment for development

python -m venv venv
source venv/bin/activate

python -m pip install -e ../bgem
python -m pip install -r requirements.txt
