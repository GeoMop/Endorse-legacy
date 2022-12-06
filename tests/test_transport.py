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
    seed = 101
    with common.workdir(f"sandbox/full_transport_{seed}", clean=False):
        # params for single container source
        source_params = dict(position=10, length=6)
        fullscale_transport(conf_file, source_params, seed)

def test_fracture_conductivity():
    common.EndorseCache.instance().expire_all()
    conf_file = os.path.join(script_dir, "test_data/config_fr_Forsmark_repo.yaml")
    cfg_fr = common.load_config(conf_file)

    rho = 998.
    g = 9.81
    visc = 0.001
    r = 10.
    K = []
    for frd in cfg_fr:
        tra = float(frd.tr_a)
        trb = float(frd.tr_b)

        T = tra * r**trb
        delta = (12*T*visc / (rho*g))**(1./3.)

        K.append(T/delta)

    print("K:",K)
    import statistics
    print("mean(K):",statistics.mean(K))

    year = 365.2425*24*3600
    # fig.5 A
    # https://onlinelibrary.wiley.com/doi/epdf/10.1111/gfl.12089
    Kf_tsx = 1 * 0.001 / year
    print("Kf_tsx", Kf_tsx)



if __name__ == "__main__":
    os.chdir(os.path.join(script_dir))
    test_macro_transport()