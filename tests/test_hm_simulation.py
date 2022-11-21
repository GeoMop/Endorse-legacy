from endorse import hm_simulation
import os
import numpy as np
from endorse import common
script_dir = os.path.dirname(os.path.realpath(__file__))


def test_run_random_samples():
    # Test execution of few random samples from the Byes inversion
    # implementation steps:
    # DONE - run single sample HM simulation equivelent to TSX model
    # - modify template to get conductivity field and porosity (possibly read results in GMSH and compute these two fields)
    # - merge config_txt and config_homogenisation
    # - run 1 sample, HM simulation but with modified geometry (2.2m diameter borehole)
    # - run 1 sample, HM simulation but with modified geometry, interpolation to other 3D tranport mesh
    # - run 1 sample, HM simulation but with modified geometry, interpolation to other 3D tranport mesh, homogenized transport execution (under 10mins localy)
    #common.EndorseCache.instance().expire_all()
    np.random.seed(0)
    conf_file = os.path.join(script_dir, "test_data/config_tsx.yaml")
    cfg = common.load_config(conf_file)
    hm_simulation.run_random_samples(cfg, 1)