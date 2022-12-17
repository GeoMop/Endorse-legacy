import os
import sys
import csv
import time
import yaml
import shutil
import aux_functions

from flow123d_simulation import endorse_2Dtest


def setup(output_dir, can_overwrite, clean):
    # create and cd workdir
    rep_dir = os.path.dirname(os.path.abspath(__file__))
    work_dir = output_dir

    # Files in the directory are used by each simulation at that level
    common_files_dir = os.path.join(work_dir, "common_files")
    # Create working directory if necessary)
    aux_functions.force_mkdir(common_files_dir, force=clean)
    os.chdir(work_dir)

    # test if config exists, copy from rep_dir if necessary
    config_file = os.path.join(common_files_dir, "config.yaml")
    if not os.path.exists(config_file):
        shutil.copyfile(os.path.join(rep_dir, "config.yaml"), config_file)

    # read config file and setup paths
    with open(config_file, "r") as f:
        config_dict = yaml.safe_load(f)

    config_dict["work_dir"] = work_dir
    config_dict["script_dir"] = rep_dir
    config_dict["common_files_dir"] = common_files_dir

    config_dict["_aux_flow_path"] = config_dict["local"]["flow_executable"].copy()
    config_dict["_aux_gmsh_path"] = config_dict["local"]["gmsh_executable"].copy()

    # copy common files
    for f in config_dict["copy_files"]:
        filepath = os.path.join(common_files_dir, f)
        if not os.path.isfile(filepath) or can_overwrite:
            shutil.copyfile(os.path.join(rep_dir, f), filepath)

    return config_dict


if __name__ == "__main__":

    # RUN THE MCMC SIMULATION
    # default parameters
    output_dir = None
    csv_data = None

    len_argv = len(sys.argv)
    assert len_argv > 2, "Specify output dir and parameters in csv file!"
    if len_argv > 1:
        output_dir = os.path.abspath(sys.argv[1])
    if len_argv > 2:
        file = sys.argv[2]
        if os.path.exists(file):
            csv_data = os.path.abspath(sys.argv[2])
        else:
            raise Exception("Missing parameters file '" + file + "'.")

    # setup paths and directories
    config_dict = setup(output_dir, can_overwrite=False, clean=False)

    print("Reading parameters from CSV: ", csv_data)
    with open(csv_data, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_NONNUMERIC)
        header = next(reader)  # skip the header
        parameters = list()
        for row in reader:
            # process each row
            parameters.append(row)

    print(parameters)

    # TODO: select random params

    # JUST RUN FLOW123D
    config_dict["solver_id"] = 0
    sim = endorse_2Dtest(config_dict)
    sim.set_parameters(data_par=parameters[0][1:])
    res, obs_data = sim.get_observations()
    print("Flow123d res: ", res)
