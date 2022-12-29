import os
import pytest
from endorse import common
from endorse.scripts.endorse_mlmc import FullScaleTransport
import subprocess
script_dir = os.path.dirname(os.path.realpath(__file__))
endorse_dir = os.path.join(script_dir, "../..")



def run_script(args):
    script_args = ['python', os.path.join(endorse_dir, 'src/endorse/scripts/endorse_mlmc.py')]
    # TODO: run as subprocess
    workdir = os.path.join(script_dir, '../sandbox/mlmc_run')

    cfg = common.load_config('inputs/config.yaml', collect_files=True)
    inputs = cfg._file_refs
    with common.workdir(workdir, inputs):
        subprocess.run(script_args + args)


# collect samples
@pytest.mark.skip
def test_FullScaleTransport_run():
    #common.EndorseCache.instance().expire_all()
    case='edz_pos02'
    #case='edz_pos10'
    #case='noedz_pos02'
    #case='noedz_pos10'

    argv = [ 'run', f'sandbox/{case}', '--clean', '--debug']
    args = FullScaleTransport.get_arguments(argv)
    pr = FullScaleTransport(f"test_data/cfg_mlmc_{case}.yaml", args)



@pytest.mark.skip
def test_script_run():
    run_script(['plot', 'cases', '*', '2 10'])

def test_script_sample():
    run_script(['run', '*', '2 10'])
