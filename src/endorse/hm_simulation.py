from typing import *
import pandas
import os
import numpy as np

import endorse.mesh_class
from .common import File, sample_from_population, workdir, dotdict
from .flow123d_simulation import endorse_2Dtest

# 2D cross-section mesh is in xy plane with center in zero
# target mesh cross-section is in yz plane
class TunnelInterpolator:
    def __init__(self, cfg_geom: dotdict, mesh: endorse.mesh_class.Mesh):
        self._mesh = mesh
        self._cfg_geom = cfg_geom
        bars = mesh.el_barycenters()
        self._barycenters = (bars[0,:], bars[1,:])

    # Maps 'target_point' from the 3D tunnel to 'point' in 2D tunnel cross-section.
    def map_point(self, target_point):
        shift = np.array([0, 0, self._cfg_geom.borehole.z_pos])
        shifted_point = target_point - shift
        # transform x<-y, y<-z
        point = (shifted_point[1], shifted_point[2])
        return point

    def interpolate_field(self, field_name, target_point, time):
        field_values = self._mesh.get_p0_values(field_name, time=time)
        point = self.map_point(target_point)
        import scipy
        # 2D cross-section mesh is in xy plane with center in zero
        val = scipy.interpolate.griddata(self._barycenters, field_values, point, method='linear')
        return val



def run_hm_simulation(config_dict: dotdict, i_sim: int, parameters: Dict[str,Union[int, float]]):
    # RUN THE MCMC SIMULATION
    # default parameters
    # output_dir = None
    # csv_data = None
    #
    # len_argv = len(sys.argv)
    # assert len_argv > 2, "Specify output dir and parameters in csv file!"
    # if len_argv > 1:
    #     output_dir = os.path.abspath(sys.argv[1])
    # if len_argv > 2:
    #     file = sys.argv[2]
    #     if os.path.exists(file):
    #         csv_data = os.path.abspath(sys.argv[2])
    #     else:
    #         raise Exception("Missing parameters file '" + file + "'.")

    # setup paths and directories
    # create and cd workdir
    work_dir_name = f"sample_{i_sim:00d}"
    #files_to_copy["common_files/config.yaml"] = "rep_dir/config.yaml"
    # TODO: seems that we need to copy the config.yaml just for reference, rather write out cfg just before simulation call ?
    rep_dir = os.path.dirname(os.path.abspath(__file__))



    work_dir = os.path.abspath(work_dir_name)
    config_dict["work_dir"] = work_dir
    config_dict["script_dir"] = rep_dir
    config_dict["common_files_dir"] = os.path.join(work_dir, "common_files")
    # TODO: make simulation independent of these three variables,
    # first two should not be necessary, the last one should be hardwired

    #config_dict["_aux_flow_path"] = config_dict["local"]["flow_executable"].copy()
    #config_dict["_aux_gmsh_path"] = config_dict["local"]["gmsh_executable"].copy()

    # copy common files
    files_to_copy = [
        (src, os.path.join("common_files", os.path.basename(src))) for src in config_dict.tsx_hm_model.copy_files
    ]

    with workdir(work_dir_name, files_to_copy):
        config_dict["solver_id"] = 0
        sim = endorse_2Dtest(config_dict)
        sim.set_parameters(parameters)
        res, obs_data = sim.get_observations()
        print("Flow123d res: ", res)


def read_bayes_sample_parameteres(parameter_file:File) -> pandas.DataFrame:
    return pandas.read_csv(parameter_file.path, dtype={'N': 'int'})


def run_random_samples(cfg, n_samples):
    df = read_bayes_sample_parameteres(File(cfg.tsx_hm_model.bayes_samples_input_file))
    i_samples = sample_from_population(n_samples, df['N'])
    for i in i_samples:
        sample_param_dict = df[i: i + 1].to_dict('records')[0]
        run_hm_simulation(cfg, i, sample_param_dict)