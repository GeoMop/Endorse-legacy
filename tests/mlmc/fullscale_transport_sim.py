import os
import copy
import shutil
import numpy as np
from typing import *

import mlmc.random.correlated_field as cf
from typing import List
from mlmc.sim.simulation import Simulation
from mlmc.quantity.quantity_spec import QuantitySpec
from mlmc.level_simulation import LevelSimulation
from endorse.fullscale_transport import compute_fields, set_source_limits, fracture_map, fullscale_transport_mesh

from endorse import common
from endorse.common import dotdict, memoize, File, call_flow, workdir, report
from endorse.mesh_class import Mesh


class FullScaleTransportSim(Simulation):

    def __init__(self, config):
        """
        :param config: Dict, simulation configuration
        """
        #super().__init__()
        self._config = config

    def level_instance(self, fine_level_params: List[float], coarse_level_params: List[float]) -> LevelSimulation:
        """
        Called from mlmc.Sampler, it creates single instance of LevelSimulation (mlmc.level_simulation)
        :param fine_level_params: fine simulation step at particular level
        :param coarse_level_params: coarse simulation step at particular level
        :return: mlmc.LevelSimulation object
        """
        config = copy.deepcopy(self._config)
        # Set sample specific parameters
        # config["fine"] = {}
        # config["coarse"] = {}
        # config["fine"]["n_steps"] = fine_level_params[0]
        # config["coarse"]["n_steps"] = coarse_level_params[0]
        # config["res_format"] = self.result_format()

        return LevelSimulation(config_dict=config,
                               calculate=FullScaleTransportSim.calculate,
                               task_size=config.transport_fullscale.mesh.fracture_mesh_step,  # @TODO: set size
                               need_sample_workspace=True)

    @staticmethod
    def calculate(cfg, seed):
        """
        Calculate fine and coarse sample and also extract their results
        :param cfg: dictionary containing simulation configuration
        :param seed: random number generator seed
        :return: np.ndarray, np.ndarray
        """
        ###################
        ### fine sample ###
        ###################
        cfg_fine = cfg.transport_fullscale
        full_mesh_file, fractures = fullscale_transport_mesh(cfg, cfg_fine.mesh, seed)

        full_mesh = Mesh.load_mesh(full_mesh_file, heal_tol=1e-4)
        el_to_fr = fracture_map(full_mesh, fractures)
        # mesh_modified_file = full_mesh.write_fields("mesh_modified.msh2")
        # mesh_modified = Mesh.load_mesh(mesh_modified_file)

        # copy files to sample dir
        shutil.copy(os.path.join(cfg["work_dir"], cfg_fine.piezo_head_input_file), os.getcwd())
        shutil.copy(os.path.join(cfg["work_dir"], cfg_fine.conc_flux_file), os.getcwd())

        input_fields_file = compute_fields(cfg, full_mesh, el_to_fr)
        large_model = File(os.path.basename(cfg_fine.piezo_head_input_file))
        conc_flux = File(os.path.basename(cfg_fine.conc_flux_file))

        params = cfg_fine.copy()
        new_params = dict(
            mesh_file=input_fields_file,
            piezo_head_input_file=large_model,
            conc_flux_file=conc_flux,
            input_fields_file=input_fields_file
        )
        params.update(new_params)
        params.update(set_source_limits(cfg))
        template = os.path.join(common.flow123d_inputs_path, cfg_fine.input_template)

        #cfg.flow_env["flow_executable"] = ["docker", "run", "-v", "{}:{}".format(os.getcwd(), os.getcwd()), "flow123d/flow123d-gnu:3.9.1", "flow123d"]
        fo = common.call_flow(cfg.flow_env, template, params)
        # set fine result
        fine_res = fo.hydro

        #####################
        ### coarse sample ###
        #####################
        coarse_res = 0

        return fine_res, coarse_res

    def result_format(self) -> List[QuantitySpec]:
        """
        Result format
        :return:
        """
        spec1 = QuantitySpec(name="conductivity", unit="m", shape=(1, 1), times=[1], locations=['0'])
        # spec2 = QuantitySpec(name="width", unit="mm", shape=(2, 1), times=[1, 2, 3], locations=['30', '40'])
        return [spec1]
