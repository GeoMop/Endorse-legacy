"""
Main module used by:
- main script: process_script.py
- fine sample script: fine_sample.py
- coarse sample script: coarse_sample.py
"""
import os
import sys
import shutil
import yaml
import numpy as np
import time as t
from typing import Any, List
import attr

src_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(src_path, '../MLMC/src'))
sys.path.append(os.path.join(src_path, '../dfn/src'))

import mlmc.sample as sample
from mlmc.simulation import Simulation
from mlmc import flow_mc
from mlmc.sample import Sample

