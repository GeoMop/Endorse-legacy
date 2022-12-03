import os
import sys
import numpy as np
from mlmc.estimator import Estimate, estimate_n_samples_for_target_variance
from mlmc.sampler import Sampler
from mlmc.sample_storage import Memory
from mlmc.sampling_pool import OneProcessPool
from mlmc.sampling_pool_pbs import SamplingPoolPBS
from mlmc.quantity.quantity import make_root_quantity
from mlmc.quantity.quantity_estimate import moments, estimate_mean
from mlmc.moments import Legendre, Monomial
from mlmc.sample_storage_hdf import SampleStorageHDF
from mlmc.plot.plots import Distribution
from mlmc.tool.process_base import ProcessBase
from mlmc import estimator
from mlmc.quantity.quantity_estimate import estimate_mean, moments
from fullscale_transport_sim import FullScaleTransportSim
from endorse import common

"""
tested parameters: run ../ --clean
"""


class FullScaleTransport:

    def __init__(self):
        args = ProcessBase.get_arguments(sys.argv[1:])

        self.work_dir = os.path.abspath(args.work_dir)
        # Add samples to existing ones
        self.clean = args.clean
        # Remove HDF5 file, start from scratch
        self.debug = args.debug
        # 'Debug' mode is on - keep sample directories
        self.use_pbs = False
        # Use PBS sampling pool
        self.n_levels = 1
        self.n_moments = 25
        self._quantile = 1e-3
        step_range = [25, 5]  # @TODO: set meaningful values or specify level parameters in a different way.
        self.level_parameters = estimator.determine_level_parameters(self.n_levels, step_range)

        # Determine number of samples at each level
        self.n_samples = estimator.determine_n_samples(self.n_levels)

        if args.command == 'run':
            self.run()
        elif args.command == 'recollect':
            self.run(recollect=True)
        elif args.command == "process":
            self.process()
        else:
            self.clean = False
            self.run(renew=True) if args.command == 'renew' else self.run()

    def run(self, renew=False, recollect=False):
        """
        Run MLMC
        :param renew: If True then rerun failed samples with same sample id
        :return: None
        """
        # Create working directory if necessary
        os.makedirs(self.work_dir, mode=0o775, exist_ok=True)

        if self.clean:
            # Remove HFD5 file
            if os.path.exists(os.path.join(self.work_dir, "mlmc_{}.hdf5".format(self.n_levels))):
                os.remove(os.path.join(self.work_dir, "mlmc_{}.hdf5".format(self.n_levels)))

        # Create sampler (mlmc.Sampler instance) - crucial class that actually schedule samples
        sampler = self.setup_config()

        # Schedule samples
        if recollect:
            raise NotImplementedError("Not supported in the released version yet")
        else:
            self.generate_jobs(sampler, n_samples=[5], renew=renew)
            self.all_collect(sampler)  # Check if all samples are finished

    def setup_config(self):
        """
        Simulation dependent configuration
        :return: mlmc.sampler instance
        """
        self.set_environment_variables()

        # Create sampling pool
        sampling_pool = self.create_sampling_pool()

        # General simulation config
        conf_file = os.path.join(self.work_dir, "test_data/config_homogenisation.yaml")
        cfg = common.load_config(conf_file)
        cfg['work_dir'] = self.work_dir
        # cfg.flow_env["flow_executable"] = ["docker", "run", "-v", "{}:{}".format(os.getcwd(), os.getcwd()),
        #                                    "flow123d/flow123d-gnu:3.9.1", "flow123d"]

        # Create simulation factory, instance of class that inherits from mlmc.sim.simulation
        simulation_factory = FullScaleTransportSim(config=cfg)

        # Create HDF sample storage
        sample_storage = SampleStorageHDF(file_path=os.path.join(self.work_dir, "mlmc_{}.hdf5".format(self.n_levels)))

        # Create sampler, it manages sample scheduling and so on
        sampler = Sampler(sample_storage=sample_storage, sampling_pool=sampling_pool, sim_factory=simulation_factory,
                          level_parameters=self.level_parameters)

        return sampler

    def set_environment_variables(self):
        """
        Set pbs config, flow123d, gmsh
        :return: None
        """
        root_dir = os.path.abspath(self.work_dir)
        while root_dir != '/':
            root_dir, tail = os.path.split(root_dir)

        if tail == 'storage' or tail == 'auto':
            # Metacentrum
            self.sample_sleep = 30
            self.init_sample_timeout = 600
            self.sample_timeout = 60
            self.adding_samples_coef = 0.1
        else:
            # Local
            self.sample_sleep = 1
            self.init_sample_timeout = 60
            self.sample_timeout = 60
            self.adding_samples_coef = 0.1

    def create_sampling_pool(self):
        """
        Initialize sampling pool, object which
        :return: None
        """
        if not self.use_pbs:
            return OneProcessPool(work_dir=self.work_dir, debug=self.debug)  # Everything runs in one process

        # Create PBS sampling pool
        sampling_pool = SamplingPoolPBS(work_dir=self.work_dir, debug=self.debug)
        # sampling_pool = OneProcessPool(work_dir=self.work_dir, debug=self.debug)

        pbs_config = dict(
            n_cores=1,
            n_nodes=1,
            select_flags={'cgroups': 'cpuacct', 'scratch_local': '10gb'},
            mem='1Gb',
            queue='charon',
            pbs_name='flow123d',
            walltime='72:00:00',
            optional_pbs_requests=[],  # e.g. ['#PBS -m ae', ...]
            home_dir='/storage/liberec3-tul/home/martin_spetlik/',
            python='python3.8',
            env_setting=['cd $MLMC_WORKDIR',
                         'module load python/3.8.0-gcc',
                         'source env/bin/activate',
                         'module use /storage/praha1/home/jan-hybs/modules',
                         'module load flow123d',
                         'module unload python-3.6.2-gcc',
                         'module unload python36-modules-gcc'],
            scratch_dir=None
        )

        sampling_pool.pbs_common_setting(flow_3=True, **pbs_config)

        return sampling_pool

    def generate_jobs(self, sampler, n_samples=None, renew=False, target_var=None):
        """
        Generate level samples
        :param n_samples: None or list, number of samples for each level
        :param renew: rerun failed samples with same random seed (= same sample id)
        :return: None
        """
        if renew:
            sampler.ask_sampling_pool_for_samples()
            sampler.renew_failed_samples()
            sampler.ask_sampling_pool_for_samples(sleep=self.sample_sleep, timeout=self.sample_timeout)
        else:
            if n_samples is not None:
                sampler.set_initial_n_samples(n_samples)
            else:
                sampler.set_initial_n_samples()
            sampler.schedule_samples()
            sampler.ask_sampling_pool_for_samples(sleep=self.sample_sleep, timeout=self.sample_timeout)
            self.all_collect(sampler)

            if target_var is not None:
                root_quantity = make_root_quantity(storage=sampler.sample_storage,
                                                   q_specs=sampler.sample_storage.load_result_format())

                moments_fn = self.set_moments(root_quantity, sampler.sample_storage, n_moments=self.n_moments)
                estimate_obj = estimator.Estimate(root_quantity, sample_storage=sampler.sample_storage,
                                                  moments_fn=moments_fn)

                # New estimation according to already finished samples
                variances, n_ops = estimate_obj.estimate_diff_vars_regression(sampler._n_scheduled_samples)
                n_estimated = estimator.estimate_n_samples_for_target_variance(target_var, variances, n_ops,
                                                                               n_levels=sampler.n_levels)

                # Loop until number of estimated samples is greater than the number of scheduled samples
                while not sampler.process_adding_samples(n_estimated, self.sample_sleep, self.adding_samples_coef,
                                                         timeout=self.sample_timeout):
                    # New estimation according to already finished samples
                    variances, n_ops = estimate_obj.estimate_diff_vars_regression(sampler._n_scheduled_samples)
                    n_estimated = estimator.estimate_n_samples_for_target_variance(target_var, variances, n_ops,
                                                                                   n_levels=sampler.n_levels)

    def set_moments(self, quantity, sample_storage, n_moments=5):
        true_domain = estimator.Estimate.estimate_domain(quantity, sample_storage, quantile=0.01)
        return Legendre(n_moments, true_domain)

    def all_collect(self, sampler):
        """
        Collect samples
        :param sampler: mlmc.Sampler object
        :return: None
        """
        running = 1
        while running > 0:
            running = 0
            running += sampler.ask_sampling_pool_for_samples(sleep=self.sample_sleep, timeout=1e-5)
            print("N running: ", running)

    def process(self):
        sample_storage = SampleStorageHDF(file_path=os.path.join(self.work_dir, "mlmc_{}.hdf5".format(self.n_levels)))
        sample_storage.chunk_size = 1024
        result_format = sample_storage.load_result_format()
        root_quantity = make_root_quantity(sample_storage, result_format)

        conductivity = root_quantity['conductivity']
        time = conductivity[1]  # times: [1]
        location = time['0']  # locations: ['0']
        values = location[0, 0]  # result shape: (1, 1)

        # Create estimator for quantities
        x_estimator = self.create_estimator(values, sample_storage)

        root_quantity_estimated_domain = Estimate.estimate_domain(root_quantity, sample_storage, self._quantile)
        root_quantity_moments_fn = Legendre(self.n_moments, root_quantity_estimated_domain)

        # There is another possible approach to calculating all moments at once and then choose quantity
        moments_quantity = moments(root_quantity, moments_fn=root_quantity_moments_fn, mom_at_bottom=True)
        moments_mean = estimate_mean(moments_quantity)
        conductivity = root_quantity['conductivity']
        time = conductivity[1]  # times: [1]
        location = time['0']  # locations: ['0']
        value_moments = location[0, 0]  # result shape: (1, 1)

        # true_domain = [-10, 10]  # keep all values on the original domain
        # central_moments_fn = Monomial(n_moments, true_domain, ref_domain=true_domain, mean=moments_mean())
        # central_moments_quantity = moments(root_quantity, moments_fn=central_moments_fn, mom_at_bottom=True)
        # central_moments_mean = estimate_mean(central_moments_quantity)
        # print("central moments mean ", central_moments_mean())

        FullScaleTransport._approx_distribution(x_estimator, self.n_levels, tol=1e-8)

    def create_estimator(self, quantity, sample_storage):
        estimated_domain = Estimate.estimate_domain(quantity, sample_storage, quantile=self._quantile)
        moments_fn = Legendre(self.n_moments, estimated_domain)
        # Create estimator for your quantity
        return Estimate(quantity=quantity, sample_storage=sample_storage,
                                              moments_fn=moments_fn)

    @staticmethod
    def _approx_distribution(estimator, n_levels, tol=1.95):
        """
        Probability density function approximation
        :param estimator: mlmc.estimator.Estimate instance, it contains quantity for which the density is approximated
        :param n_levels: int, number of MLMC levels
        :param tol: Tolerance of the fitting problem, with account for variances in moments.
        :return: None
        """
        distr_obj, result, _, _ = estimator.construct_density(tol=tol)
        distr_plot = Distribution(title="distributions", error_plot=None)
        distr_plot.add_distribution(distr_obj)

        if n_levels == 1:
            samples = estimator.get_level_samples(level_id=0)[..., 0]
            distr_plot.add_raw_samples(np.squeeze(samples))
        distr_plot.show(None)
        distr_plot.reset()

    @staticmethod
    def determine_level_parameters(n_levels, step_range):
        """
        Determine level parameters,
        In this case, a step of fine simulation at each level
        :param n_levels: number of MLMC levels
        :param step_range: simulation step range
        :return: list of lists
        """
        assert step_range[0] > step_range[1]
        level_parameters = []
        for i_level in range(n_levels):
            if n_levels == 1:
                level_param = 1
            else:
                level_param = i_level / (n_levels - 1)
            level_parameters.append([step_range[0] ** (1 - level_param) * step_range[1] ** level_param])
        return level_parameters


if __name__ == "__main__":
    pr = FullScaleTransport()

