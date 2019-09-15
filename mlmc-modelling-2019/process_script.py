"""
Main MLMC execution stript.
"""

import sys
import os
import yaml
import shutil
import numpy as np
import matplotlib.pyplot as plt
src_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(src_path, '../MLMC/src'))
sys.path.append(os.path.join(src_path, '../gmsh-tools/src'))

from mlmc import mlmc
from mlmc import base_process
from mlmc.estimate import Estimate, CompareLevels
from mlmc.moments import Legendre, Monomial

import pbs
import fracture_model as fm
from mlmc.simulation import Simulation
from mlmc import flow_mc
from mlmc.sample import Sample

def load_config_dict(cfg):
    with open(cfg, "r") as f:
        return yaml.safe_load(f)



class SimulationSample(Sample):
    """
    Object produced by the Simulation object.
    - Used to compute both fine and coarse samples.
    - Has to be serializable

    """
    def __init__(self, sample_dir, sample_id, seed, fine):
        """
        :param config: Configuration dict of the sample:
        seed:
        sample_dir:
        fine_sample:
        """

        super().__init__(
            sample_id=sample_id,
            directory=sample_dir
        )
        self.random_seed = seed
        self.fine = fine

    def queue_execution(self, pbs):
        sample_config = dict(
            random_seed=self.random_seed,
            level=
            h_step)
        if not self.fine:
            return
        self.job_id = package_dir


    def extract_result(self):
        result_file = ["FINISHED_COARSE.npy", "FINISHED_FINE.npy"][self.fine]
        return np.load(result_file)









class FractureFlowSimulation(Simulation):
    total_sim_id = 0

    def __init__(self, mesh_step, level_id=None, config=None, clean=False, parent_fine_sim=None):
        if level_id is not None:
            self.sim_id = level_id
        else:
            self.sim_id = FractureFlowSimulation.total_sim_id
            FractureFlowSimulation.total_sim_id += 1

        self.step = mesh_step
        # Pbs script creater
        self.pbs_creater = config['pbs']
        self.n_fine_elements = 0


        # Prepare base workdir for this mesh_step
        output_dir = config['output_dir']
        self.work_dir = os.path.join(output_dir, 'sim_%d_step_%f' % (self.sim_id, self.step))
        flow_mc.force_mkdir(self.work_dir, clean)


        self.coarse_sim = None
        self.coarse_sim_set = False

    @property
    def is_fine_sim(self):
        return self.coarse_sim_set


    def n_ops_estimate(self):
        """
        Number of operations
        :return: int
        """
        return (1000 / self.step)**2

    def set_coarse_sim(self, coarse_sim=None):
        """
        Set coarse simulation ot the fine simulation so that the fine can generate the
        correlated input data sample for both.

        Here in particular set_points to the field generator
        :param coarse_sim
        """
        self.coarse_sim = coarse_sim
        self.coarse_sim_set = True

    def generate_random_sample(self):
        """
        Generate common sample input for fine and coare simulation.
        Here nothing, all is done in fine_sim.simulation_sample
        """

    def write_sample_config(self, sample_dir):
        if self.coarse_sim:
            coarse_step = self.coarse_sim.step
        else:
            coarse_step = 0.0
        sample_config = dict(
            seed=np.random.randint(0, np.iinfo(np.int32).max),
            do_coarse=self.coarse_sim is not None,
            h_fine_step=self.step,
            h_coarse_step=coarse_step,
            config_path=os.path.join(src_path, "config.yaml")
        )
        with open(os.path.join(sample_dir, "sample_config.yaml"), "w") as f:
            yaml.dump(sample_config, f)

    def simulation_sample(self, sample_tag, sample_id, start_time=0):
        """

        :param sample_tag:
        :param sample_id:
        :param start_time:
        :return:
        """


        if not self.is_fine_sim:
            sample_dir = ""
            package_dir=""
        else:
            # Fine sim.
            sample_dir = os.path.join(self.work_dir, sample_tag)
            flow_mc.force_mkdir(sample_dir)
            for f in ['flow_templ.yaml']:
                shutil.copy(os.path.join(src_path, f), os.path.join(sample_dir, f))
            self.write_sample_config(sample_dir)
            # Fine sample starts execution job for both samples
            lines = [
                'cd {sample_dir}',
                'python3 {src_dir}/both_sample.py sample_config.yaml 2>&1 | tee both_sample_out',
            ]
            package_dir = self.pbs_creater.add_realization(
                weight=self.n_fine_elements,
                lines=lines,
                sample_dir=sample_dir,
                src_dir=src_path)

        sample = Sample(directory=sample_dir, sample_id=sample_id,
                      job_id=package_dir)
        sample.is_fine_sample = self.is_fine_sim
        return sample



    def _extract_result(self, sample):
        if sample.is_coarse_sample:
            pass
        else:
            pass



class Process(base_process.Process):

    def create_pbs_object(self, output_dir, clean):
        """
        Initialize object for PBS execution
        :param output_dir: Output directory
        :param clean: bool, if True remove existing files
        :return: None
        """
        pbs_work_dir = os.path.join(output_dir, "scripts")
        num_jobs = 0
        if os.path.isdir(pbs_work_dir):
            num_jobs = len([_ for _ in os.listdir(pbs_work_dir)])

        self.pbs_obj = pbs.Pbs(pbs_work_dir,
                               job_count=num_jobs,
                               qsub=self.pbs_config['qsub'],
                               clean=clean)
        self.pbs_obj.pbs_common_setting(**self.pbs_config)

    def setup_config(self, n_levels, clean):
        """
        Simulation dependent configuration
        :param n_levels: Number of levels
        :param clean: bool, If True remove existing files
        :return: mlmc.MLMC instance
        """
        self.step_range = (100, 10)
        # Set pbs config, flow123d, gmsh, ...
        self.config_dict = load_config_dict(os.path.join(src_path, "config.yaml"))

        root_dir = os.path.abspath(self.work_dir)
        while root_dir != '/':
            root_dir, tail = os.path.split(root_dir)

        self.flow123d = self.config_dict["flow_executable"]
        self.pbs_config = dict(
            job_weight=250000,  # max number of elements per job
            n_cores=1,
            n_nodes=1,
            select_flags=[],
            mem='8gb',
            queue='charon',
            qsub=None)
        if(self.config_dict["metacentrum"]):
            self.sample_sleep = 30
            self.init_sample_timeout = 600
            self.sample_timeout = 0
            self.pbs_config.update(dict(
                qsub='/usr/bin/qsub',
                modules=[

                ]
            ))
        else:
            # pbs_config is necessary for the local run but is not used
            # as long as the pbs key is set to None
            self.pbs_config['qsub'] = None
            self.sample_sleep = 1
            self.init_sample_timeout = 60
            self.sample_timeout = 60

        self.mc_samples = self.config_dict["mc_samples"]
        output_dir = os.path.join(self.work_dir, "output_{}".format(n_levels))

        reuse_samples = self.config_dict.get('reuse_samples', None)
        if reuse_samples:
            clean = False
            self.config_dict['metacentrum'] = False
            self.config_dict['finish_sleep'] = 0

        # remove existing files
        if clean:
            self.rm_files(output_dir)

        # copy files work_dir
        if src_path != self.work_dir:
            for f in ["config.yaml", "flow_templ.yaml"]:
                orig = os.path.join(src_path, f)
                work = os.path.join(self.work_dir, f)
                if clean or not os.path.isfile(work):
                    shutil.copy(orig, work)

        # Init pbs object
        self.create_pbs_object(output_dir, clean)
        # if (self.config_dict["metacentrum"]):
        #
        # else:
        #     self.pbs_obj = None

        simulation_config = dict(
            pbs=self.pbs_obj,
            output_dir=output_dir
            )




        self.options['output_dir'] = output_dir
        mlmc_obj = mlmc.MLMC(n_levels, FractureFlowSimulation.factory(self.step_range, config=simulation_config, clean=clean),
                                  self.step_range, self.options)

        if clean or reuse_samples:
            # Create new execution of mlmc
            mlmc_obj.create_new_execution()
        else:
            # Use existing mlmc HDF file
            mlmc_obj.load_from_file()
        return mlmc_obj


    def run(self):
        """
        Run mlmc
        :return: None
        """
        os.makedirs(self.work_dir, mode=0o775, exist_ok=True)
        self.n_moments = 10

        mlmc_list = []
        # Run one level Monte-Carlo method
        for nl in [2]:  # , 2, 3, 4,5, 7, 9]:
            mlmc = self.setup_config(nl, clean=True)

            # self.n_sample_estimate(mlmc)
            self.generate_jobs(mlmc, n_samples=[self.mc_samples])
            mlmc_list.append(mlmc)

        self.all_collect(mlmc_list)



    def process(self):
        """
        Use collected data
        :return: None
        """
        assert os.path.isdir(self.work_dir)
        mlmc_est_list = []

        for nl in [1]:  # high resolution fields
            mlmc = self.setup_config(nl, clean=False)
            # Use wrapper object for working with collected data
            mlmc_est = Estimate(mlmc)
            mlmc_est_list.append(mlmc_est)

        print("PROCESS FINISHED :)")




if __name__ == "__main__":
    process = Process()
