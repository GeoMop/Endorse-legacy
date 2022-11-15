import os.path
from typing import *
import yaml
import subprocess
import shutil
from pathlib import Path
import numpy as np




class workdir:
    """
    Context manager for creation and usage of a workspace dir.

    name: the workspace directory
    inputs: list of files and directories to copy into the workspaceand
        TODO: fine a sort of robust ad portable reference
    clean: if true the workspace would be deleted at the end of the context manager.
    """
    def __init__(self, name:str="sandbox", inputs:List[str] = None, clean=False):
        if inputs is None:
            inputs = []
        self._inputs = inputs
        self.work_dir = os.path.abspath(name)
        Path(self.work_dir).mkdir(parents=True, exist_ok=True)
        self._clean = clean
        self._orig_dir = os.getcwd()

    def __enter__(self):
        for f in self._inputs:
            if os.path.isdir(f):
                shutil.copytree(f, self.work_dir, dirs_exist_ok=True)
            else:
                shutil.copy2(f, self.work_dir)
        os.chdir(self.work_dir)

        return self.work_dir

    def __exit__(self, type, value, traceback):
        os.chdir(self._orig_dir)
        if self._clean:
            shutil.rmtree(self.work_dir)


class dotdict(dict):
    """
    dot.notation access to dictionary attributes
    """
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

    def __getattr__(self, item):
        return self[item]

    @classmethod
    def create(cls, cfg : Any):
        """
        - recursively replace all dicts by the dotdict.
        """
        if isinstance(cfg, dict):
            items = ( (k, cls.create(v)) for k,v in cfg.items())
            return dotdict(items)
        elif isinstance(cfg, list):
            return [cls.create(i) for i in cfg]
        elif isinstance(cfg, tuple):
            return tuple([cls.create(i) for i in cfg])
        else:
            return cfg


def load_config(path):
    """
    Load configuration from given file replace, dictionaries by dotdict
    """
    with open(path) as f:
        cfg = yaml.safe_load(f)
        return dotdict.create(cfg)


def substitute_placeholders(file_in: str, file_out: str, params: Dict[str, Any]):
    """
    In the template `file_in` substitute the placeholders of format '<name>'
    according to the dict `params`. Write the result to `file_out`.
    """
    used_params = []
    with open(file_in, 'r') as src:
        text = src.read()
    for name, value in params.items():
        placeholder = '<%s>' % name
        n_repl = text.count(placeholder)
        if n_repl > 0:
            used_params.append(name)
            text = text.replace(placeholder, str(value))
    with open(file_out, 'w') as dst:
        dst.write(text)
    return used_params


def check_conv_reasons(log_fname):
    """
    Check correct convergence of the solver.
    Reports the divergence reason and returns false in case of divergence.
    """
    with open(log_fname, "r") as f:
        for line in f:
            tokens = line.split(" ")
            try:
                i = tokens.index('convergence')
                if tokens[i + 1] == 'reason':
                    value = tokens[i + 2].rstrip(",")
                    conv_reason = int(value)
                    if conv_reason < 0:
                        print("Failed to converge: ", conv_reason)
                        return False
            except ValueError:
                continue
    return True

# Directory for all flow123d main input templates.
# These are considered part of the software.
_script_dir = os.path.dirname(os.path.realpath(__file__))
flow123d_inputs_path = os.path.join(_script_dir, "../flow123d_inputs")

def call_flow(cfg, file_in, params):
    """
    Run Flow123d in actual work dir with main input given be given template and dictionary of parameters.

    1. prepare the main input file from filebase_in + "_tmpl.yamlL"
    2. run Flow123d
    """
    in_dir, template = os.path.split(file_in)
    suffix = "_tmpl.yaml"
    assert template[-len(suffix):] == suffix
    filebase = template[:-len(suffix)]
    main_input = filebase + ".yaml"
    substitute_placeholders(file_in, main_input, params)

    arguments = cfg.flow_executable.copy()
    arguments.append(main_input)
    print("Running: ", " ".join(arguments))
    with open(filebase + "_stdout", "w") as stdout:
        with open(filebase + "_stderr", "w") as stderr:
            completed = subprocess.run(arguments, stdout=stdout, stderr=stderr)
            # print(completed)
    print("Exit status: ", completed.returncode)
    success = completed.returncode == 0
    if not success:
        with open(filebase + "_stderr", "r") as stderr:
            print(stderr.read())
        raise Exception("Flow123d ended with error")
    conv_check = check_conv_reasons(os.path.join("output", "flow123.0.log"))
    print("converged: ", conv_check)
    return success  # and conv_check




# TODO: running with stdout/ stderr capture, test for errors, log but only pass to the main in the case of
# true error


def sample_from_population(n_samples:int, frequency:Union[np.array, int]):
    if type(frequency) is int:
        frequency = np.full(len(frequency), 1, dtype=int)
    else:
        frequency = np.array(frequency, dtype=int)

    cumul_freq = np.cumsum(frequency)
    total_samples = np.sum(frequency)
    samples = np.random.randint(0, total_samples, size=n_samples + 1)
    samples[-1] = total_samples # stopper
    sample_seq = np.sort(samples)
    # put samples into bins given by cumul_freq
    bin_samples = np.empty_like(samples)
    i_sample = 0
    for ifreq, c_freq in enumerate(cumul_freq):

        while sample_seq[i_sample] < c_freq:
            bin_samples[i_sample] = ifreq
            i_sample += 1

    return bin_samples[:-1]