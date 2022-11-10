from endorse import common
import os

script_dir = script_dir = os.path.dirname(os.path.realpath(__file__))


def test_load_config():
    conf_file = os.path.join(script_dir, "test_data/config.yaml")
    cfg = common.load_config(conf_file)
    assert isinstance(cfg.geometry, common.dotdict)
    assert isinstance(cfg.geometry, common.dotdict)