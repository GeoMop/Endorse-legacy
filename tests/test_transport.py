import pytest
import os

#import endorse.macro_flow_model
from endorse import common
from endorse.fullscale_transport import fullscale_transport, input_files


script_dir = os.path.dirname(os.path.realpath(__file__))

#@pytest.mark.skip
def test_macro_transport():
   # with common.workdir("sandbox"):
    common.EndorseCache.instance().expire_all()
    conf_file = os.path.join(script_dir, "test_data/config_homo_tsx.yaml")
    #cfg = common.load_config(conf_file)
    #files = input_files(cfg.transport_fullscale)
    seed = 19
    with common.workdir(f"sandbox/full_transport_{seed}", clean=False):
        # params for single container source
        source_params = dict(position=2, length=4)
        fullscale_transport(conf_file, source_params, seed)

if __name__ == "__main__":
    os.chdir(os.path.join(script_dir))
    test_macro_transport()