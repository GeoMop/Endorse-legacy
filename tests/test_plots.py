import os
import pytest
from endorse import plots, common
from bgem.gmsh.gmsh_io import GmshIO
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

script_dir = os.path.dirname(os.path.realpath(__file__))


@pytest.mark.skip
def test_plot_mesh():
    mesh = GmshIO('test_data/rectangle_2x5.msh')
    fig, ax = plt.subplots(figsize=(10,10))
    plots.plot_mesh_2d(ax, mesh)
    X = np.linspace(0,1,100)
    Y = X**2
    #ax.plot(X,Y)
    plt.show()

    print('DONE')

def test_plot_source():
    conf_file = os.path.join(script_dir, "test_data/config_homo_tsx.yaml")
    cfg = common.load_config(conf_file)
    cfg_fine = cfg.transport_fullscale
    conc_flux = common.File(cfg_fine.conc_flux_file)
    plots.plot_source(conc_flux)


@pytest.mark.skip
def test_plot_errorbar():
    data_dict = {"endz_10": "path to edz_pos10/mlmc_1.hdf5",
                 "noedz_10": "path to noedz_pos10/mlmc_1.hdf5",
                 "endz_02": "path to  edz_pos02/mlmc_1.hdf5",
                 "noedz_02": "path to noedz_pos02/mlmc_1.hdf5"
                 }
    plots.plot_quantile_errorbar(data_dict)
