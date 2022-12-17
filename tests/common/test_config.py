from endorse import common
import os
import numpy as np
from scipy import stats
import pytest

script_dir = script_dir = os.path.dirname(os.path.realpath(__file__))




def test_dotdict():
    pass

def test_apply_variant():
    cfg = dict(
        a = [0, 1],
        b = dict(
            a='a',
            b='b'
        )
    )
    cfg = common.dotdict(cfg)
    variant = {
        'a/1' : 2,
        'b/b' : 'c'
    }
    new_cfg = common.apply_variant(cfg, variant)
    assert new_cfg['a'][0] == 0
    assert new_cfg['a'][1] == 2
    assert new_cfg['b']['a'] == 'a'
    assert new_cfg['b']['b'] == 'c'

    # test errors

    # missing key in dict
    with pytest.raises(KeyError):
        common.apply_variant(cfg, {'c':0})
    # integer key for dict
    with pytest.raises(KeyError):
        common.apply_variant(cfg, {'0':0})
    # empty path
    with pytest.raises(KeyError):
        common.apply_variant(cfg, {'':0})
    # missing index of list
    with pytest.raises(IndexError):
        common.apply_variant(cfg, {'a/2':0})
    # indexing list by key
    with pytest.raises(IndexError):
        common.apply_variant(cfg, {'a/b':0})

def test_load_config():
    conf_file = os.path.join(script_dir, "../test_data/config.yaml")
    cfg = common.load_config(conf_file)
    assert isinstance(cfg.geometry, common.dotdict)
    assert isinstance(cfg.geometry, common.dotdict)