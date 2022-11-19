from typing import *
import logging
import os
import attrs
from . import memoize, File, report, substitute_placeholders, workdir
import subprocess

def search_file(basename, extensions):
    """
    Return first found file or None.
    """
    if type(extensions) is str:
        extensions = (extensions,)
    for ext in extensions:
        if os.path.isfile(basename + ext):
            return File(basename + ext)
    return None

@attrs.define
class EquationOutput:
    @staticmethod
    def from_cwd(basenames):
        spatial_base, balance_base, observe_base = basenames
        return EquationOutput(
            spatial=search_file(spatial_base+"_fields", (".msh", ".pvd")),
            balance=search_file(balance_base+"_balance", ".txt"),
            observe=search_file(observe_base+"_observe", ".txt")
        )
    spatial: File
    balance: File
    observe: File


class FlowOutput:

    def __init__(self, process: subprocess.CompletedProcess, stdout, stderr):
        self.process = process
        self.stdout = File(stdout.name)
        self.stderr = File(stderr.name)
        with workdir("output"):
            self.log = File("flow123.0.log")
            self.hydro = EquationOutput.from_cwd(("flow", "flow", "flow"))
            self.solute = EquationOutput.from_cwd(("solute", "solute", "solute"))
            self.mechanic = EquationOutput.from_cwd(("mechanics", "mechanics", "mechanics"))


    def check_conv_reasons(self):
        """
        Check correct convergence of the solver.
        Reports the divergence reason and returns false in case of divergence.
        """
        with open(self.log.path, "r") as f:
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


@report
@memoize
def call_flow(cfg:'dotdict', file_in:File, params: Dict[str,str]) -> FlowOutput:
    """
    Run Flow123d in actual work dir with main input given be given template and dictionary of parameters.

    1. prepare the main input file from filebase_in + "_tmpl.yamlL"
    2. run Flow123d

    TODO: pass only flow configuration
    """
    in_dir, template = os.path.split(file_in)
    suffix = "_tmpl.yaml"
    assert template[-len(suffix):] == suffix
    filebase = template[:-len(suffix)]
    main_input = filebase + ".yaml"
    substitute_placeholders(file_in, main_input, params)

    arguments = cfg.flow_executable.copy()
    arguments.append(main_input)
    logging.info("Running Flow123d: " + " ".join(arguments))
    with open(filebase + "_stdout", "w") as stdout:
        with open(filebase + "_stderr", "w") as stderr:
            completed = subprocess.run(arguments, stdout=stdout, stderr=stderr)
    fo = FlowOutput(completed, stdout, stderr)

    logging.info(f"Exit status: {completed.returncode}")
    if completed.returncode != 0:
        with open(filebase + "_stderr", "r") as stderr:
            print(stderr.read())
        raise Exception("Flow123d ended with error")
    conv_check = fo.check_conv_reasons()
    logging.info(f"converged: {conv_check}")
    return fo

# TODO:
# - call_flow variant with creating dir, copy,


