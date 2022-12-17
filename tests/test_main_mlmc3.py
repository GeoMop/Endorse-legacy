import pytest
from endorse import common
from endorse.mlmc.mlmc_main import FullScaleTransport

# collect samples
def test_FullScaleTransport_run():
    #common.EndorseCache.instance().expire_all()

    #case='edz_pos02'
    #case='edz_pos10'
    #case='noedz_pos02'
    case='noedz_pos10'

    argv = [ 'run', f'sandbox/{case}', '--clean', '--debug']
    args = FullScaleTransport.get_arguments(argv)
    pr = FullScaleTransport(f"test_data/cfg_mlmc_{case}.yaml", args)


def test_script():
    # TODO: run as subprocess
    pass