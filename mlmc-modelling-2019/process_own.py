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

import pbs


"""
Script overview:
- independent of mlmc
- schedule unfinished samples
- wait for unfinished samples
- move failed (should copy them to avoid their reschdulling
- collect results, compute statistics (when all samples are finished, checked in extract_results)
"""








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
            i_level=self.i_level,
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
            print("Schedule: ", sample_dir)
            os.makedirs(sample_dir, mode=0o775, exist_ok=True)
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
                print("FAILED -------------------------------------------------\n",
                      sample_dir)
                traceback.print_exc()
                print("-----------------------------")
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
            self.
            #self.test_homogenity()
            self.test_homogenity_xy()
            self.test_isotropy()
            #self.test_isotropy_alt()
            self.compute_variogram()
            self.calculate_field_parameters()

    @property
    def cond_xy(self):
        if self._cond_xy is None:
            self._cond_xy = np.concatenate(self.cond_field_xy, axis=0)
        return self._cond_xy

    @property
    def cond_tn(self):
        if self._cond_tn is None:
            self._cond_tn = np.concatenate(self.cond_field_values, axis=0)
        return self._cond_tn

    def test_homogenity(self):
        points = self.cond_xy
        cond_tn = self.cond_tn
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
        points = self.cond_xy
        cond_tn = self.cond_tn
        mean_c_tn = np.mean(cond_tn, axis=0)
        print(mean_c_tn)
        fig, axes = plt.subplots(nrows=2, ncols=3)

        for iy, axes_x in enumerate(axes):
            for ix, ax in enumerate(axes_x):
                coord = [(0,0), (1,1), (0,1)][ix]
                X = points[:, iy]
                Y = cond_tn[:, coord[0], coord[1]]
                sc = ax.scatter(X, Y, s=1)
                Y0 = np.full_like(X, mean_c_tn[coord[0], coord[1]])
                ax.plot(X, Y0, c='red')
                ax.set_ylim([np.min(Y), np.max(Y)])
                if iy == 0:
                    ax.set_title(["Cxx", "Cyy", "Cxy"][ix])
                if ix == 0:
                    ax.set_ylabel("avg. over " + ["Y", "X"][iy])
        fig.suptitle("Homogenity of the conductivity field")
        fig.savefig("homogenity.pdf")
        plt.close(fig)

    def test_var_homogenity(self):
        points = self.cond_xy
        cond_tn = self.cond_tn
        mean_c_tn = np.mean(cond_tn, axis=0)
        print(mean_c_tn)
        fig, axes = plt.subplots(nrows=2, ncols=3)

        for iy, axes_x in enumerate(axes):
            for ix, ax in enumerate(axes_x):
                coord = [(0,0), (1,1), (0,1)][ix]
                X = points[:, iy]
                Y = cond_tn[:, coord[0], coord[1]]
                sc = ax.scatter(X, Y, s=1)
                Y0 = np.full_like(X, mean_c_tn[coord[0], coord[1]])
                ax.plot(X, Y0, c='red')
                ax.set_ylim([np.min(Y), np.max(Y)])
                if iy == 0:
                    ax.set_title(["Cxx", "Cyy", "Cxy"][ix])
                if ix == 0:
                    ax.set_ylabel("avg. over " + ["Y", "X"][iy])
        fig.suptitle("Homogenity of the conductivity field")
        fig.savefig("homogenity.pdf")
        plt.close(fig)


    def test_isotropy(self):
        cond_tn = self.cond_tn
        mean_c_tn = np.mean(cond_tn, axis=0)
        print("Mean tensor: \n", mean_c_tn)
        e_val, e_vec = np.linalg.eigh(mean_c_tn)
        print("Cmin: {} Cmax: {} angle: {}"
            .format(e_val[0], e_val[1], np.angle(e_vec[0, 1] + e_vec[1, 1] * 1j)))

        fig, ax = plt.subplots(nrows=1, ncols=1)
        angle = []
        c_val = []
        for ctn in cond_tn:
            e_val, e_vec = np.linalg.eigh(ctn)
            angle.append(np.angle(e_vec[0, :] + e_vec[1, :] * 1j))
            c_val.append(e_val)
            angle.append(np.angle(-e_vec[0, :] - e_vec[1, :] * 1j))
            c_val.append(e_val)
        angle = np.array(angle)
        c_val = np.array(c_val)
        #ax.set_aspect('equal')
        ax.scatter(angle[:, 0]/np.pi, c_val[:, 0], label='C min eigv', s=0.5, c='blue')
        ax.scatter(angle[:, 1]/np.pi, c_val[:, 1],  label='C max eigv', s=0.5, c='darkblue')

        angle = np.linspace(-np.pi, np.pi, 200)
        vec = np.stack([np.cos(angle), np.sin(angle)], axis=0)
        mean_flux = mean_c_tn @ vec

        ax.plot(angle/np.pi, np.linalg.norm(mean_flux, axis=0),
                label='C mean', color='orange', linewidth=2)
        Y = np.full_like(angle, 1e-10)
        ax.plot(angle/np.pi, Y,
                label='C bulk', color='red', linewidth=2)
        ax.xaxis.set_major_formatter(plt.FormatStrFormatter('%g $\pi$'))
        import matplotlib
        ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(base=0.25))
        ax.set_ylim(np.min(c_val.flatten()), np.max(c_val.flatten()))
        ax.set_yscale('log')
        ax.legend()
        fig.suptitle("Isotropy  of the conductivity field")
        fig.savefig("isotropy.pdf")
        plt.close(fig)

    def test_isotropy_alt(self):
        x_size = 100
        y_size = 100
        cond_tn = self.cond_tn
        angle = np.linspace(0, 2*np.pi, x_size)
        vec = np.stack([np.cos(angle), np.sin(angle)], axis=1)
        cond = []
        for v in vec:
            c = cond_tn[np.random.choice(len(cond_tn), y_size), :, :]
            cond.append(np.log10(np.linalg.norm(np.dot(v, c), axis=1)))
        cond = np.concatenate(cond)
        fig, ax = plt.subplots(nrows=1, ncols=1)
        angle = np.repeat(angle, y_size)
        ax.scatter(angle.flatten(), cond.flatten(), s=0.5)
        self.set_x_pi_ticks(ax)
        ax.set_ylabel("|C @ unit_vector|")
        ax.set_xlabel("Unit vector angle")
        ax.legend()
        fig.suptitle("Isotropy  of the conductivity field")
        fig.savefig("isotropy_alt.pdf")
        plt.close(fig)

    def set_x_pi_ticks(self, ax):
        from matplotlib import ticker
        ax.xaxis.set_major_formatter(ticker.FuncFormatter(
            lambda x, pos: '{:.2f} $\pi$'.format(x/np.pi)

        ))
        ax.xaxis.set_major_locator(ticker.MultipleLocator(base=0.25*np.pi))


    def compute_variogram(self):
        point_samples = self.cond_xy
        cond_samples = self.cond_tn
        max_length = 1000

        # Select pairs to sample various point distances
        #radius = 0.5 * np.linalg.norm(points.max_pt - points.min_pt, 2)
        n_samples = len(point_samples)
        assert len(cond_samples) == n_samples
        n_pairs_per_sample = 10000
        dist_list = []
        variogram_list = []
        for points, conds in zip(point_samples, cond_samples):
            var = np.var(conds, axis=0, ddof=1)
            n_vals = len(points)
            assert len(conds) == n_vals
            pairs = np.random.choice(n_vals, (n_pairs_per_sample, 2))
            pair_dists = np.linalg.norm(points[pairs[:, 0]] - points[pairs[:, 1]], axis=1)
            pair_variogram = 0.5 * np.abs(conds[pairs[:, 0]] - conds[pairs[:, 1]]) ** 2 / var
            dist_list.append(pair_dists)
            variogram_list.append(pair_variogram)
        dists = np.concatenate(dist_list)
        variograms = np.concatenate(variogram_list)
        indices = np.argsort(dists)
        dists = dists[indices]
        n_reliable_dist = np.argmin(dists < max_length/2)
        dists = dists[:n_reliable_dist]
        variograms = variograms[indices[:n_reliable_dist]]

        n_cells = 10
        breaks = np.linspace(0, len(dists), n_cells, endpoint=False, dtype=int)


        breaks = list(breaks) + [len(dists)]
        print("cell sizes: \n",[breaks[i+1] - breaks[i]
                      for i in range(len(breaks) - 1)])
        print("dist intervals:\n", [dists[breaks[i+1]-1] - dists[breaks[i]]
                      for i in range(len(breaks) - 1)])
        cell_dists = [np.mean(dists[breaks[i]: breaks[i + 1]])
                      for i in range(len(breaks) - 1)]
        cell_variogram = [np.mean(variograms[breaks[i]: breaks[i+1]], axis=0)
                          for i in range(len(breaks) - 1)]
        cell_variogram = np.array(cell_variogram)

        fig, axes = plt.subplots(nrows=2, ncols=2)
        for iy, axes_x in enumerate(axes):
            for ix, ax in enumerate(axes_x):
                #ax.scatter(dists, variograms[:, iy, ix], s=0.5)
                label = "C{}{}".format(["x", "y"][ix], ["x", "y"][iy])
                ax.plot(cell_dists, cell_variogram[:, iy, ix],
                        color='red', label=label)
                ax.set_xlabel("corr. length")
                ax.legend()
        fig.suptitle("Empirical variogram.")
        fig.savefig("correlations.pdf")
        plt.close(fig)



    def calculate_field_parameters(self):
        """
        Assuming:
        - uncorrelated tensor field
        - tensor given by Cmax, Cmin, angle
        - we assume that [log(Cmax), log(Cmin)] is Bivariete normal distr
        - angle is uniform on [0, pi]
        - Cpair is independent on angle
        We compute Cmax, Cmin, angle samples (! for Cmax~Cmin the angle is not defined)
        so we drop such samples for angle related statistics
        :return:
        """
        cond_tn = self.cond_tn
        print("n samples: ", len(cond_tn))
        angle = []  # angle of larger eig val
        c_val = []  # pairs of eigvals [cmin, cmax]
        for ctn in cond_tn:
            e_val, e_vec = np.linalg.eigh(ctn)
            a = np.angle(e_vec[0, 1] + e_vec[1, 1] * 1j)
            if a < 0:
                a += np.pi
            angle.append(a)
            c_val.append(np.log10(e_val))
            #c_val.append(e_val)
        angle = np.array(angle)
        c_val = np.array(c_val)

        # abs conductivity condition
        #i_reliable_angle = c_val[:, 1] / c_val[:, 0] > 2
        # log conductivity condition
        i_reliable_angle = c_val[:, 1] - c_val[:, 0] > np.log10(2)
        print("reliable samples: ", sum(i_reliable_angle))

        # reliable angle ECDF plot
        rel_angle = angle[i_reliable_angle]
        # Cval - angle correlation
        rel_c_val = c_val[i_reliable_angle, :]
        rel_c_angle = np.concatenate((rel_c_val, rel_angle[:,None]), axis=1)
        c_angle_cov = np.cov(rel_c_angle.T)
        # print("cmin, cmax, angle cov. matrix:\n", c_angle_cov)
        # print("cmin, cmax, angle corr. matrix:\n",
        #       c_angle_cov /
        #       np.sqrt((c_angle_cov.diagonal()[:, None] * c_angle_cov.diagonal())[None, :])
        #       )
        c_angle_corr = np.corrcoef(rel_c_angle.T)
        print("cmin, cmax, angle corr. matrix:\n", c_angle_corr)

        fig, ax = plt.subplots(nrows=1, ncols=1)
        ax.plot(np.sort(rel_angle), np.linspace(0, 1, len(rel_angle)))
        self.set_x_pi_ticks(ax)
        ax.set_ylabel("empirical CDF")
        ax.set_xlabel("maximum eigen vector angle $\phi$")
        #fig.suptitle()

        ax.text(0, 0.9,
                "n. samples: {},\n"
                "with eigenvalue ratio > 2".format(len(rel_angle)))
        ax.text(0.4*np.pi, 0,
                """
                Correlations:
                $\log(C_{{max}}), \phi: \quad\quad\quad$ {:.3f}
                $\log(C_{{min}}), \phi: \quad\quad\quad$ {:.3f}
                $\log(C_{{max}}), \log(C_{{min}}):$ {:.3f}
                """.format(c_angle_corr[1, 2],c_angle_corr[0, 2],c_angle_corr[0, 1])
                )
        fig.savefig("angle_ecdf.pdf")
        plt.close(fig)

        c_val = rel_c_val
        mean_cval = np.mean(c_val, axis=0)
        print("Mean eigen values:\n", mean_cval)
        cov_cval = np.cov(c_val.T)
        print("Covarriance eigen values:\n", cov_cval)

        fig, axes = plt.subplots(nrows=1, ncols=2)

        cov_eval, cov_egvec = np.linalg.eigh(cov_cval)
        inv_cov_sqrt = (cov_egvec @ np.diag(1.0 / np.sqrt(cov_eval))).T

        std_cval = inv_cov_sqrt @ (c_val.T - mean_cval[:, None])
        from scipy.stats import norm
        for iax, ax in enumerate(axes):
            Y = norm.ppf(np.linspace(0, 1, len(c_val)))
            ax.plot(np.sort(std_cval[iax, :]), Y)
            ax.set_xlabel("samples")
            ax.set_ylabel("Q log_norm")
        fig.savefig("QQ_conductivity.pdf")
        
    def eigenvals_correlation(self):
        pass


    def calculate_field_params_mcmc(self):
        import pymc3
        cov_sqrt = (cov_egvec @ np.diag(np.sqrt(cov_eval))).T
        with pymc3.Model() as model:
            cond_indep = pymc3.Normal("conv_indep", mu=0, shape=2)
            cond_dep = cov_sqrt @ (c_indep) + log_cond_mean
            unrotated_tn = np.diag(cond_dep)
            angle = pymc3.Uniform("conv_angle", lower=0, upper=2 * np.pi)
            c, s = np.cos(angle), np.sin(angle)
            rot_mat = np.array([[c, -s], [s, c]])
            cond_2d = rot_mat.T @ unrotated_tn @ rot_mat



class Process():
    def __init__(self, work_dir):
        self.work_dir = os.path.abspath(work_dir)
        os.makedirs(self.work_dir, mode=0o775, exist_ok=True)
        work_config_path = os.path.join(self.work_dir, "config.yaml")
        if not os.path.exists(work_config_path):
            shutil.copy(os.path.join(src_path, "config.yaml"), self.work_dir)
        with open(work_config_path, "r") as f:
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
        if self.config_dict['metacentrum']:
            pbs_config['qsub'] = 'qsub'
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
            sim_step = level_config['step']
            sim = FractureFlowSimulation(
                il, self.pbs, sim_step,
                level_config['n_samples'], last_step, self.work_dir)
            self.levels.append(sim)
            last_step = sim_step
                               
        for l in reversed(self.levels):
            l.run_level()

    def move_failed(self, failed):
        failed_dir = os.path.join(self.work_dir, "FAILED")
        os.makedirs(failed_dir, mode=0o775, exist_ok=True)
        for sample_dir in failed:
            try:
                shutil.move(sample_dir, failed_dir)
            except shutil.Error:
                pass
            # make empty finished dir
            os.makedirs(sample_dir, mode=0o775, exist_ok=True)
            with open(os.path.join(sample_dir, "FINISHED"), "w") as f:
                f.write('done')


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
    work_dir = sys.argv[1]
    process = Process(work_dir)
    process.run()
    process.wait()
