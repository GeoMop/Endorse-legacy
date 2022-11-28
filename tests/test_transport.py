import pytest
import os

import endorse.macro_flow_model
from endorse import common
from endorse.fullscale_transport import fullscale_transport, input_files


script_dir = os.path.dirname(os.path.realpath(__file__))

#@pytest.mark.skip
def test_macro_transport():
   # with common.workdir("sandbox"):
    common.EndorseCache.instance().expire_all()
    conf_file = os.path.join(script_dir, "test_data/config_homogenisation.yaml")
    cfg = common.load_config(conf_file)
    files = input_files(cfg.transport_fullscale)
    seed = 100
    with common.workdir("sandbox/full_transport", inputs=files, clean=False):
        # params for single container source
        source_params = dict(position=10, length=6)
        fullscale_transport(cfg, source_params, seed)
