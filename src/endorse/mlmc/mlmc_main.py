import os
import sys
import numpy as np
from mlmc.estimator import Estimate
from mlmc.sampler import Sampler
from mlmc.sampling_pool import OneProcessPool, ProcessPool
from mlmc.sampling_pool_pbs import SamplingPoolPBS
from mlmc.quantity.quantity import make_root_quantity
from mlmc.moments import Legendre
from mlmc.sample_storage_hdf import SampleStorageHDF
from mlmc.plot.plots import Distribution
from mlmc import estimator
from mlmc.quantity.quantity_estimate import estimate_mean, moments
from mlmc.sim.fullscale_transport_sim import FullScaleTransportSim

from endorse.common.config import load_config
import os

_script_dir = os.path.dirname(os.path.realpath(__file__))
_endorse_repository = os.path.abspath(os.path.join(_script_dir, '../../../'))


"""
tested parameters: run ../ --clean
"""


class FullScaleTransport:

    def __init__(self, main_cfg_file, args):

        self.work_dir = os.path.abspath(args.work_dir)
        self.cfg_file = main_cfg_file
        cfg = load_config(os.path.join(main_cfg_file))
        self.cfg = cfg

        #cfg.flow_env.mlmc.singularity
        #self.singularity_image = os.path.abspath(cfg.mlmc.singularity_image)
        #self.endorse_repository = os.path.abspath(args.endorse_dir)
        # Add samples to existing ones
        self.clean = args.clean
        # Remove HDF5 file, start from scratch
        self.debug = args.debug
        # 'Debug' mode is on - keep sample directories
        self.use_pbs = True
        # Use PBS sampling pool
        self.n_levels = 1
        self.n_moments = 25
        self._quantile = 1e-3
        step_range = [25, 5]  # @TODO: set meaningful values or specify level parameters in a different way.
        self.level_parameters = estimator.determine_level_parameters(self.n_levels, step_range)

        # Determine number of samples at each level
        self.n_samples = [ 20 ]

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
            self.generate_jobs(sampler, n_samples=self.n_samples, renew=renew)
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
        # conf_file = os.path.join(self.work_dir, "test_data/config_homogenisation.yaml")
        # cfg = self.load_config(conf_file)
        config = {}
        config['work_dir'] = self.work_dir
        config["flow_executable"] = ["flow123d"]
        config['source_params'] = dict(position=10, length=6)
        config['main_cfg_file'] = os.path.abspath(self.cfg_file)
        config["mesh_steps"] = {self.level_parameters[0][0]: 50} # @TODO: check values
        #config = dict(work_dir=self.work_dir)
        # cfg.flow_env["flow_executable"] = ["docker", "run", "-v", "{}:{}".format(os.getcwd(), os.getcwd()),
        #                                    "flow123d/flow123d-gnu:3.9.1", "flow123d"]
        #cfg.flow_env["flow_executable"] = ["flow123d"]


        # Create simulation factory, instance of class that inherits from mlmc.sim.simulation
        simulation_factory = FullScaleTransportSim(config=config)

        # Create HDF sample storage
        sample_storage = SampleStorageHDF(file_path=os.path.join(self.work_dir, "mlmc_{}.hdf5".format(self.n_levels)))

        # Create sampler, it manages sample scheduling and so on
        sampler = Sampler(sample_storage=sample_storage, sampling_pool=sampling_pool, sim_factory=simulation_factory,
                          level_parameters=self.level_parameters)

        return sampler

    # def load_config(self, path):
    #     """
    #     Load configuration from given file replace, dictionaries by dotdict
    #     uses pyyaml-tags namely for:
    #     include tag:
    #         geometry: <% include(path="config_geometry.yaml")>
    #     """
    #     #YamlIncludeConstructor.add_to_loader_class(loader_class=yaml.FullLoader, base_dir=os.path.dirname(path))
    #     with open(path) as f:
    #         cfg = yaml.load(f, Loader=yaml.FullLoader)
    #     print("cfg ", cfg)
    #     return cfg#dotdict.create(cfg)

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
            #self.init_sample_timeout = 600
            self.sample_timeout = 60
            self.adding_samples_coef = 0.1
        else:
            # Local
            self.sample_sleep = 1
            #self.init_sample_timeout = 60
            self.sample_timeout = 60
            self.adding_samples_coef = 0.1

    def create_sampling_pool(self):
        """
        Initialize sampling pool, object which
        :return: None
        """
        cfg_pbs = self.cfg.mlmc.get('pbs', None)
        if cfg_pbs is None:
            return OneProcessPool(work_dir=self.work_dir, debug=self.debug)  # Everything runs in one process
            #return ProcessPool(4, work_dir=self.work_dir, debug=self.debug)  # Everything runs in one process

        # Create PBS sampling pool
        sampling_pool = SamplingPoolPBS(work_dir=self.work_dir, debug=self.debug)

        singularity_img = self.cfg.mlmc.singularity_image
        pbs_config = dict(
            optional_pbs_requests=[],  # e.g. ['#PBS -m ae', ...]
            #home_dir="Why we need the home dir!! Should not be necessary.",
            #home_dir='/storage/liberec3-tul/home/martin_spetlik/',
            python=f'singularity exec {singularity_img} venv/bin/python3',
            #python='singularity exec {} /usr/bin/python3'.format(self.singularity_path),
            env_setting=[#'cd $MLMC_WORKDIR',
                         "export SINGULARITY_TMPDIR=$SCRATCHDIR",
                         "export PIP_IGNORE_INSTALLED=0",
                         'cd {}'.format(_endorse_repository),
                         'singularity exec {} ./setup.sh'.format(singularity_img),
                         'singularity exec {} venv/bin/python3 -m pip install scikit-learn'.format(singularity_img)
                         #'module load python/3.8.0-gcc',
                         #'source env/bin/activate',
                         #'module use /storage/praha1/home/jan-hybs/modules',
                         #'module load flow123d',
                         #'module unload python-3.6.2-gcc',
                         #'module unload python36-modules-gcc'
                         ],
            scratch_dir=None
        )
        pbs_config.update(cfg_pbs)
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

        conductivity = root_quantity['indicator_conc']
        time = conductivity[1]  # times: [1]
        location = time['0']  # locations: ['0']
        values = location[0, 0]  # result shape: (10, 1)

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


# class dotdict(dict):
#     """
#     dot.notation access to dictionary attributes
#     TODO: keep somehow reference to the original YAML in order to report better
#     KeyError origin.
#     """
#     __setattr__ = dict.__setitem__
#     __delattr__ = dict.__delitem__
#
#     def __getattr__(self, item):
#         try:
#             return self[item]
#         except KeyError:
#             return self.__getattribute__(item)
#
#     @classmethod
#     def create(cls, cfg : Any):
#         """
#         - recursively replace all dicts by the dotdict.
#         """
#         if isinstance(cfg, dict):
#             items = ( (k, cls.create(v)) for k,v in cfg.items())
#             return dotdict(items)
#         elif isinstance(cfg, list):
#             return [cls.create(i) for i in cfg]
#         elif isinstance(cfg, tuple):
#             return tuple([cls.create(i) for i in cfg])
#         else:
#             return cfg
    @staticmethod
    def get_arguments(arguments):
        """
        Getting arguments from console
        :param arguments: list of arguments
        :return: namespace
        """
        import argparse
        parser = argparse.ArgumentParser()

        parser.add_argument('command', choices=['run', 'collect', 'renew'],
                            help='run - create new execution,'
                                 'collect - keep collected, append existing HDF file'
                                 'renew - renew failed samples, run new samples with failed sample ids (which determine random seed)')
        parser.add_argument('work_dir', help='Work directory')
        #parser.add_argument('singularity_path', help='Path to singularity image')
        #parser.add_argument('endorse_dir', help='Path to endorse repository')
        parser.add_argument("-c", "--clean", default=False, action='store_true',
                            help="Clean before run, used only with 'run' command")
        parser.add_argument("-d", "--debug", default=False, action='store_true',
                            help="Keep sample directories")

        args = parser.parse_args(arguments)

        return args

if __name__ == "__main__":
    args = FullScaleTransport.get_arguments(sys.argv[1:])
    pr = FullScaleTransport("config_homo_tsx.yaml", args)

