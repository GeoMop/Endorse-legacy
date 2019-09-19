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
#sys.path.append(os.path.join(src_path, '../gmsh-tools/src'))

from mlmc import mlmc
from mlmc import base_process
from mlmc.estimate import Estimate, CompareLevels
from mlmc.moments import Legendre, Monomial

from mlmc import pbs
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

        self.cond_field_xy = []
        self.cond_field_values = []
        self.running_samples = {}
        self.finished_samples = {}
        # Register produced samples. After all are collected we perform averaging.

        self.result_struct = [["value", "cxx", "cxy", "cyy"],
                              ["f8", "f8", "f8", "f8"]]


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

    def simulation_sample(self, i_samplesample_tag, sample_id, start_time=0):
        """

        :param sample_tag:
        :param sample_id:
        :param start_time:
        :return:
        """


        if not self.is_fine_sim:
            sample = Sample(
                        directory=self._last_fine_sample.directory,
                        sample_id=sample_id,
                        job_id=self._last_fine_sample.job_id)
            sample._fine_sample = self._last_fine_sample
        else:
            # Fine sim.
            sample_dir = os.path.join(self.work_dir, sample_tag)
            if not os.path.exists(os.path.join(sample_dir, "FINISHED")):
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
            else:
                package_dir = "finished_job"

            sample = Sample(directory=sample_dir, sample_id=sample_id,
                          job_id=package_dir)
            self.running_samples[sample_id] = sample
            if self.coarse_sim:
                self.coarse_sim._last_fine_sample = sample


        return sample




    def _extract_result(self, sample):
        finished_file = os.path.join(sample.directory, "FINISHED")
        finished = False
        if os.path.exists(finished_file):
            with open(finished_file, "r") as f:
                content = f.read().split()
            finished = len(content) == 1 and content[0] == "done"
        if finished:
            with open(os.path.join(sample.directory, "summary.yaml"), "r") as f:
                summary_dict = yaml.load(f, Loader=yaml.FullLoader)
            if self.is_fine_sim:
                cond_tn = np.array(summary_dict['fine']['cond_tn'][0])
                self.cond_field_xy.append(np.array(summary_dict['coarse_ref']['pos']))
                self.cond_field_values.append(np.array(summary_dict['coarse_ref']['cond_tn']))
            else:
                cond_tn = np.array(summary_dict['coarse']['cond_tn'][0])
            return [0, cond_tn[0, 0], cond_tn[1, 0], cond_tn[1, 1]]
        else:
            return None



    def compute_cond_field_properties(self):
        if self.cond_field_xy:
            #self.test_homogenity()
            self.test_homogenity_xy
            #self.compute_variogram(point_samples, cond_samples)

    def test_homogenity(self):
        points = np.concatenate(self.cond_field_xy, axis=0)
        cond_tn = np.concatenate(self.cond_field_values, axis=0)
        mean_c_tn = np.mean(cond_tn, axis=0)
        cond_diff = cond_tn - mean_c_tn
        print(mean_c_tn)

        fig, axes = plt.subplots(nrows=2, ncols=2)

        X, Y = points.T
        for iy, axes_x in enumerate(axes):
            for ix, ax in enumerate(axes_x):
                sc = ax.scatter(X, Y, s=1, c=cond_diff[:, ix, iy])

        fig.subplots_adjust(right=0.8)
        cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
        fig.colorbar(sc, cax=cbar_ax)

        plt.show()

    def test_homogenity_xy(self):
        points = np.concatenate(self.cond_field_xy, axis=0)
        cond_tn = np.concatenate(self.cond_field_values, axis=0)
        mean_c_tn = np.mean(cond_tn, axis=0)
        cond_diff = cond_tn - mean_c_tn
        print(mean_c_tn)

        fig, axes = plt.subplots(nrows=2, ncols=3)

        X, Y = points.T
        for iy, axes_x in enumerate(axes):
            for ix, ax in enumerate(axes_x):
                coord = [(0,0), (0,1), (1,1)][ix]
                sc = ax.scatter(points[:, iy], cond_tn[:, coord[0], coord[1]], s=1)
                ax.plot(points[:, iy], mean_c_tn[coord[0], coord[1]], c='red')


        fig.subplots_adjust(right=0.8)
        cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
        fig.colorbar(sc, cax=cbar_ax)

        plt.show()

    def compute_variogram(self, point_samples, cond_samples):

        # Select pairs to sample various point distances
        #radius = 0.5 * np.linalg.norm(points.max_pt - points.min_pt, 2)
        n_samples = len(point_samples)
        assert len(cond_samples) == n_samples
        n_pairs_per_sample = 100
        dist_list = []
        variogram_list = []
        for points, conds in zip(point_samples, cond_samples):
            n_vals = len(points)
            assert len(conds) == n_vals
            pairs = np.random.choice(n_vals, (n_pairs_per_sample, 2))

            pair_dists = np.linalg.norm(points.points[pairs[:, 0]] - points.points[pairs[:, 1]], axis=1)
            pair_variogram = np.abs(conds[pairs[:, 0]] - conds[pairs[:, 1]]) ** 2

            indices = np.argsort(pair_dists)
            pair_dists = pair_dists[indices]
            pair_variogram = pair_variogram[indices]
            dist_list.append(pair_dists)
            variogram_list.append(pair_variogram)
        dists = np.concatenate(dist_list)
        variograms = np.concatenate(variogram_list)
        n_cells = 10
        breaks = np.linspace(0, len(dists), n_cells, endpoint=False, dtype=int)
        cell_dists = dists[breaks]
        cell_variogram = []
        breaks = list(breaks) + [len(dists)]
        start = 0
        for end in breaks[1:]:
            cell_variogram.append(np.mean(variograms[start, end], axis=0))
            start = end

    def plot_variogram(self, breaks, variogram_tn):
        fig = plt.figure()
        # setting the axis limits in [left, bottom, width, height]
        #rect = [0.1, 0.1, 0.8, 0.8]
        ax_xx = fig.add_subplot(2, 2, 0)
        ax_xy = fig.add_subplot(2, 2, 1)
        ax_yx = fig.add_subplot(2, 2, 2)
        ax_yy = fig.add_subplot(2, 2, 3)

        breaks = [(breaks[i] + breaks[i+1])/2 for i in range(len(breaks))]
        ax_xx.plot(breaks, variogram_tn[:, 0, 0])
        ax_xy.plot(breaks, variogram_tn[:, 1, 0])
        ax_yx.plot(breaks, variogram_tn[:, 0, 1])
        ax_yy.plot(breaks, variogram_tn[:, 1, 1])









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

    def setup(self):
        os.makedirs(self.work_dir, mode=0o775, exist_ok=True)
        levels = [
            (FractureFlowSimulation(100, level_id=0), None),
            (FractureFlowSimulation(10, level_id=1), FractureFlowSimulation(100, level_id=1))
        ]
        n_samples = [1, 3]
        for ns, lev in zip(n_samples, reversed(levels))
            fine_sim, coarse_sim = lev
            fine_sim.set_coarse_sim(coarse_sim)
            for i_sample in range(ns):
                fine_sim.sample_simulation(i_sample)

    def run(self):
        """
        Run mlmc
        :return: None
        """

        # self.n_moments = 10
        #
        # mlmc_list = []
        # # Run one level Monte-Carlo method
        # for nl in [2]:  # , 2, 3, 4,5, 7, 9]:
        #     mlmc = self.setup_config(nl, clean=False)
        #
        #     # self.n_sample_estimate(mlmc)
        #     #n_samples = nl * [self.mc_samples]
        #     self.generate_jobs(mlmc, n_samples=[1,3])
        #     mlmc_list.append(mlmc)
        #
        # self.all_collect(mlmc_list)



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
