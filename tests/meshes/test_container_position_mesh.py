# Test src/container_position_mesh.py
# TODO: rename as the container position need not to be related to the container spacing
# fine EDZ (0.3): about  300k elements, 10 s


import os
import numpy as np
from endorse import common
from endorse import container_position_mesh
#from bgem.stochastic.fracture import Fracture


script_dir = script_dir = os.path.dirname(os.path.realpath(__file__))

def test_make_mesh():
    # about 280 k elements
    # conf_file = os.path.join(script_dir, "./config_full_coarse.yaml")
    conf_file = os.path.join(script_dir, "./config_full_edz_fine.yaml")
    cfg = common.load_config(conf_file)
    #fractures = [
    #    Fracture(4, np.array([]), np.array(), )
    #]
    fractures = []
    container_position_mesh.make_mesh(cfg.geometry, fractures, 0, "test_container_position_mesh.msh")
