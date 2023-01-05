import pytest
import os

#import endorse.macro_flow_model
from endorse import common, flow123d_inputs_path
from endorse.fullscale_transport import transport_run, transport_2d


script_dir = os.path.dirname(os.path.realpath(__file__))

def test_flow123d_templates():
    template = flow123d_inputs_path.joinpath("transport_fullscale_tmpl.yaml")
    assert os.path.isfile(template)

#@pytest.mark.skip
def test_macro_transport():
   # with common.workdir("sandbox"):
    #common.EndorseCache.instance().expire_all()
    conf_file = os.path.join(script_dir, "test_data/config.yaml")
    #cfg = common.load_config(conf_file)
    #files = input_files(cfg.transport_fullscale)
    seed = 19
    with common.workdir(f"sandbox/full_transport_{seed}", clean=False):
        # params for single container source
        cfg = common.load_config(conf_file)
        cfg['transport_fullscale']['end_time'] = 120
        ind_time_max = transport_run(cfg, seed)
        print("Result:", ind_time_max)

@pytest.mark.skip
def test_transport_2d():
   # with common.workdir("sandbox"):
    common.EndorseCache.instance().expire_all()
    conf_file = os.path.join(script_dir, "test_data/config.yaml")
    seed = 19
    with common.workdir(f"sandbox/transport_2d_{seed}", clean=False):
        # params for single container source
        cfg = common.load_config(conf_file)
        cfg['transport_fullscale']['end_time'] = 120
        ind_time_max = transport_2d(cfg, seed)
        print("Result:", ind_time_max)

if __name__ == "__main__":
    os.chdir(os.path.join(script_dir))
    test_macro_transport()