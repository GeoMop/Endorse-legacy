"""
Own process script without MLMC.
"""

import sys
import os
import time
import yaml
import shutil
import traceback
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









class FractureFlowSimulation():
    total_sim_id = 0

    def __init__(self, i_level, pbs_obj, mesh_step, n_samples, coarse_step, work_dir):
        self.i_level = i_level
        self.step = mesh_step
        self.coarse_step = coarse_step
        # Pbs script creater
        self.pbs_creater = pbs_obj
        self.n_samples = n_samples

        self.cond_field_xy = []
        self.cond_field_values = []
        self.running_samples = {}
        self.finished_samples = {}
        self.work_dir = work_dir
        # Register produced samples. After all are collected we perform averaging.


    def run_level(self):
        """
        :return: [ (level, i_sample, dir) ]
        """
        level_dir = "sim_{}_step_{:.6f}".format(self.i_level, self.step)
        level_dir = os.path.join(self.work_dir, level_dir)
        os.makedirs(level_dir, mode=0o775, exist_ok=True)
        for i_sample in range(self.n_samples):
            sample_dir = "L{:02d}_F_S{:07}".format(self.i_level, i_sample)
            sample_dir = os.path.join(level_dir, sample_dir)
            os.makedirs(sample_dir, mode=0o775, exist_ok=True)
            self.simulation_sample(i_sample, sample_dir)
        return self.running_samples



    def write_sample_config(self, sample_dir):
        sample_config = dict(
            seed=np.random.randint(0, np.iinfo(np.int32).max),
            do_coarse=self.coarse_step is not None,
            h_fine_step=self.step,
            h_coarse_step=self.coarse_step,
            config_path=os.path.join(self.work_dir, "config.yaml")
        )
        config_path = os.path.join(sample_dir, "sample_config.yaml")
        if not os.path.exists(config_path):
            with open(config_path, "w") as f:
                yaml.dump(sample_config, f)

    def simulation_sample(self, i_sample, sample_dir):
        """
        :param sample_tag:
        :param sample_id:
        :param start_time:
        :return:
        """

        # Fine sim.
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
                weight=100,
                lines=lines,
                sample_dir=sample_dir,
                src_dir=src_path)
        else:
            package_dir = "finished_job"

        self.running_samples[i_sample] = sample_dir



    def extract_results(self):
        new_running = {}
        failed = []
        for i_sample, sample_dir in self.running_samples.items():
            try:
                result = self.extract_result(sample_dir)
                if result is None:
                    new_running[i_sample] = sample_dir
                    continue
                self.finished_samples[i_sample] = (sample_dir, result)
            except:
                traceback.print_exc()
                failed.append(sample_dir)
        self.running_samples = new_running
        if len(self.running_samples) == 0:
            self.compute_cond_field_properties()
        return failed

    def extract_result(self, sample_dir):
        """
        Return:
         None - not yet finished
         pair of (fine, coarse) sample, or None
        :param sample_dir:
        :return:
        """
        finished_file = os.path.join(sample_dir, "FINISHED")
        finished = False
        if os.path.exists(finished_file):
            with open(finished_file, "r") as f:
                content = f.read().split()
            finished = len(content) == 1 and content[0] == "done"
        if finished:
            with open(os.path.join(sample_dir, "summary.yaml"), "r") as f:
                summary_dict = yaml.load(f) #, Loader=yaml.FullLoader

            fine_cond_tn = np.array(summary_dict['fine']['cond_tn'][0])
            if self.coarse_step is not None:
                coarse_cond_tn = np.array(summary_dict['coarse']['cond_tn'][0])
                self.cond_field_xy.append(np.array(summary_dict['coarse_ref']['pos']))
                self.cond_field_values.append(np.array(summary_dict['coarse_ref']['cond_tn']))
            else:
                coarse_cond_tn = None
            return fine_cond_tn, coarse_cond_tn
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









class Process():
    def __init__(self, work_dir):
        self.work_dir = os.path.abspath(work_dir)
        os.makedirs(self.work_dir, mode=0o775, exist_ok=True)
        shutil.copy(os.path.join(src_path, "config.yaml"), self.work_dir)
        with open(os.path.join(self.work_dir, "config.yaml"), "r") as f:
            self.config_dict = yaml.load(f) #, Loader=yaml.FullLoader
        




    def make_pbs(self):
        pbs_config = dict(
            job_weight=250000,  # max number of elements per job
            n_cores=1,
            n_nodes=1,
            select_flags=[],
            mem='8gb',
            queue='charon',
            qsub=None)
        pbs_obj = pbs.Pbs(self.work_dir,
                               job_count=0,
                               qsub=pbs_config['qsub']
                               )
        pbs_obj.pbs_common_setting(**pbs_config)
        return  pbs_obj


    def run(self):
        os.makedirs(self.work_dir, mode=0o775, exist_ok=True)
        self.pbs=self.make_pbs()

        self.levels = []
        last_step = None
        for il, level_config in enumerate(self.config_dict['levels']):
            sim = FractureFlowSimulation(il, self.pbs, level_config['step'], level_config['n_samples'], last_step, self.work_dir)
            self.levels.append(sim)
                               
        for l in reversed(self.levels):
            l.run_level()

    def move_failed(self, failed):
        failed_dir = os.path.join(self.work_dir, "FAILED")
        os.makedirs(failed_dir, mode=0o775, exist_ok=True)
        for sample_dir in failed:
            shutil.move(sample_dir, failed_dir)


    def wait(self):
        self.pbs.execute()
        n_running = 1
        while (n_running):
            n_running = 0
            for sim in self.levels:
                failed = sim.extract_results()
                self.move_failed(failed)
                n_running += len(sim.running_samples)
            time.sleep(1)


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
    np.random.seed(123)
    process = Process("output_2")
    process.run()
    process.wait()
