import pytest
from endorse.common.config import dotdict
from endorse.mlmc.mlmc_main import FullScaleTransport

def test_FullScaleTransport():
     argv = [ 'run', 'sandbox', '--clean', '--debug']
     args = FullScaleTransport.get_arguments(argv)
     pr = FullScaleTransport("test_data/config_homo_tsx.yaml", args)


def test_script():
    # TODO: run as subprocess
    pass