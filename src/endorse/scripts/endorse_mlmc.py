from typing import *
import os
import sys
import numpy as np
import attrs
import argparse
import fnmatch
import shutil

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

from endorse import common
from endorse import plots

_script_dir = os.path.dirname(os.path.realpath(__file__))
_endorse_repository = os.path.abspath(os.path.join(_script_dir, '../../../'))

MAIN_CONFIG_FILE = 'config.yaml'

"""
tested parameters: run ../ --clean
"""


class FullScaleTransport:

    def __init__(self, main_cfg_file, args):

        self.work_dir = os.path.abspath(args.work_dir)
        self.cfg_file = main_cfg_file
        cfg = common.load_config(os.path.join(main_cfg_file))
        self.cfg = cfg

        #cfg.flow_env.mlmc.singularity
        #self.singularity_image = os.path.abspath(cfg.mlmc.singularity_image)
        #self.endorse_repository = os.path.abspath(args.endorse_dir)
        # Add samples to existing ones
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


    def run(self, renew=False, recollect=False):
        """0
        Run MLMC
        :param renew: If True then rerun failed samples with same sample id
        :return: None
        """
        # Create working directory if necessary
        os.makedirs(self.work_dir, mode=0o775, exist_ok=True)

        # Create sampler (mlmc.Sampler instance) - crucial class that actually schedule samples
        sampler = self.setup_config()

        self.generate_jobs(sampler, n_samples=self.n_samples, renew=renew)
        self.all_collect(sampler)  # Check if all samples are finished

    def recollect(self):
        self.run(recollect=True)

    def renew(self):
        self.clean = False
        self.run(renew=True)





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

CaseName = NewType('CaseName', str)
CasePatch = Dict[common.config.Path, common.dotdict]

@attrs.define
class SourceDensity:
    """
    For odd length N the density is (1/N, ..., 1/N),  N items
    for even length N the density is (1/2N, 1/N, ..., 1/N, 1/2N) N + 1 items
    """
    # first container indexed from 0
    center: int
    # number of containers in the source
    length: int

    def plot_label(self):
        return f"pos: {self.center}, len: {self.length}"

    def fs_name_items(self):
        c_str =  f"{self.center:03d}"
        if self.length > 1:
            return [c_str, str(self.source.length)]
        else:
            return [c_str]

    #@staticmethod
    #def from_center(center: int, length: int):
    #    return SourceDensity(center - , length)

@attrs.define
class SimCase:
    case_name: str
    case_patch: CasePatch
    source: SourceDensity

    @property
    def directory(self):
        items = [self.case_name, *self.source.fs_name_items()]
        return "-".join(items)

    @property
    def hdf5_path(self):
        return os.path.join(self.directory, "mlmc_1.hdf5")


    def mean_std_log(self):
        i_quantile = 1
        sample_storage = SampleStorageHDF(file_path=self.hdf5_path)
        sample_storage.chunk_size = 1024
        result_format = sample_storage.load_result_format()
        root_quantity = make_root_quantity(sample_storage, result_format)

        conductivity = root_quantity['indicator_conc']
        time = conductivity[1]  # times: [1]
        location = time['0']  # locations: ['0']
        values = location  # result shape: (10, 1)
        values = values[i_quantile, 0]  # selected quantile
        values = values.select(values < 1e-1)
        values = np.log(values)
        samples = self._get_samples(values, sample_storage)[0]


        q_mean = estimate_mean(values)
        val_squares = estimate_mean(np.power(values - q_mean.mean, 2))
        std = np.sqrt(val_squares.mean)

        return q_mean.mean[0], std[0], samples

    def _get_samples(self, quantity, sample_storage):
        n_moments = 2
        estimated_domain = Estimate.estimate_domain(quantity, sample_storage, quantile=0.001)
        moments_fn = Legendre(n_moments, estimated_domain)
        estimator = Estimate(quantity=quantity, sample_storage=sample_storage, moments_fn=moments_fn)
        samples = estimator.get_level_samples(level_id=0)[..., 0]
        return samples


@attrs.define
class SimCases:
    """
    Organize set of stochastic simulation variants.
    That consists of cartesian product of variants (cases) of forward simulation parameters
    and distribution of contamination source (sources), currently just as discrete containers.

    Set of all named configuration cases is defined by 'cases.yaml' (file must be referenced in main 'config.yaml') that
    contains dictionary {<case_name> : <dict of config changes>}.

    Source distribution should be considered as prescribed set of concentration sources, i.s. concentration source density C(x)
    is sum of densities for individual containers. The diffusion through the bentonite would be described by a time changing
    and random parameter however the distribution is independent on the container position.
    This way the process is linear with respect to the concentration density and better describes differences between
    compared configurations.

    To this end we consider source given by its position and span (how many source containers we consider) the concentration is normalized
    so that total initial concentration is 1. Time decrease is part of configuration (see cases).
    <i> = [i:i+1:1]
    <i>-<n> = [i-n/2 : i+n/2+(1) :1]
    """


    cfg_cases : Dict[CaseName, CasePatch]
    source_densities: List[SourceDensity]

    @staticmethod
    def def_args(parser):
        help = """
               Space separated names of cases. Subset of cases defined as keys of `cases.yaml`."
               """
        parser.add_argument("cases", help=help)
        help = '''
        Defines basis functions of the linear space of source densities.Space separated index ranges. E.g. "1:10:2 2 6:8"
        '''
        parser.add_argument("sources", help=help)

    @classmethod
    def initialize(cls, args):
        """
        Initialize from cfg._cases file and arguments.
        :param args:
        :return:
        """
        cfg = common.load_config(MAIN_CONFIG_FILE)
        cases = cls.active_cases(args.cases, cfg.cases_file)
        sources = cls.source_basis(args.sources, cfg.geometry.containers.n_containers)
        return cls(cases, sources)

    @staticmethod
    def active_cases(arg_cases, cases_file):
        cfg = common.load_config(cases_file)
        if not isinstance(cfg, dict):
            cfg = {}
        all_cases = list(cfg.keys())
        selected_cases = []
        for case_pattern in arg_cases:
            selected_cases.extend(fnmatch.filter(all_cases, case_pattern))
        return {k:cfg[k] for k in selected_cases}

    @staticmethod
    def source_basis(arg_sources, n_containers):
        sources = []
        all_containers = list(range(n_containers))
        for slice_token in arg_sources.split(" "):
            slice_items = [int(x.strip()) if x.strip() else None for x in slice_token.split(':')]
            if len(slice_items) == 1:
                # single number interpreted as single index contrary to the standard slice
                source_slice = slice(slice_items[0], slice_items[0] + 1, 1)
            else:
                source_slice = slice(*slice_items)
            subset = all_containers[source_slice]
            if not subset:
                continue
            step = subset[1] - subset[0] if len(subset) > 1 else 1
            for i in subset:
                sources.append(SourceDensity(i, step))
        return sources

    def iterate(self):
        for case_key, case_patch in self.cfg_cases.items():
            for source in self.source_densities:
                yield SimCase(case_key, case_patch, source)




class CleanCmd:
    @staticmethod
    def def_args(parser):
        help="Remove whole cases directories."
        parser.add_argument('--all', action='store_true', help=help)
        SimCases.def_args(parser)

    # Remove HFD5 file
    def execute(self, args):
        cases = SimCases.initialize(args)
        for case in cases.iterate():
            try:
                if args.all:
                    shutil.rmtree(case.directory, ignore_errors=False, onerror=None)
                else:
                    os.remove(os.path.join(case.directory, f"mlmc_{cases.cfg.mlmc.n_levels}.hdf5"))
            except FileNotFoundError:
                pass





@attrs.define
class RunCmd:
    @staticmethod
    def def_args(parser):
        parser.add_argument("-d", "--debug", default=False, action='store_true',
                        help="Keep sample directories")
        SimCases.def_args(parser)

    def execute(self, args):
        pass



class CasesPlot:

    @staticmethod
    def def_args(parser):
        SimCases.def_args(parser)

    def execute(self, args):
        cases = SimCases.initialize(args)
        data = [(case.case_name, case.source.plot_label(), *case.mean_std_log()) for case in cases.iterate()]
        print(data)
        plots.plot_log_errorbar_groups(data, 'conc ' + r'$[g/m^3]$')


@attrs.define
class PlotCmd:
    @staticmethod
    def def_args(parser):
        add_subparsers(parser, 'plot', 'plot_class', [CasesPlot])

    def execute(self, args):
        plot_instance = args.plot_class()
        plot_instance.execute(args)



def add_subparsers(parser:argparse.ArgumentParser, suffix, fn_arg:str, classes: List[Any]):
    subparsers = parser.add_subparsers()
    for cls in classes:
        cmd = cls.__name__.lower().rstrip(suffix)
        subparser = subparsers.add_parser(cmd)
        cls.def_args(subparser)
        subparser.set_defaults(**{fn_arg: cls})

def get_arguments(arguments):
    """
    Getting arguments from console
    :param arguments: list of arguments
    :return: namespace
    """
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-w', '--workdir',
                        default=os.getcwd(),
                        type=str, help='Main directory of the whole project. Default is current directory.')
    add_subparsers(parser, 'cmd', 'cmd_class', [CleanCmd, RunCmd, PlotCmd])
    args = parser.parse_args(arguments)
    return args


if __name__ == "__main__":
    args = get_arguments(sys.argv[1:])
    with common.workdir(args.workdir):
        command_instance = args.cmd_class()
        command_instance.execute(args)


