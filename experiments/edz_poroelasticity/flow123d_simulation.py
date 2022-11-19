import os
import subprocess

import numpy as np
import shutil
import yaml
import traceback

import matplotlib.pyplot as plt

import aux_functions


def generate_time_axis(config_dict):
    end_time = float(config_dict["end_time"])
    output_times = config_dict["output_times"]

    # create time axis
    times = []
    for dt in output_times:
        b = float(dt["begin"])
        s = float(dt["step"])
        e = float(dt["end"])
        times.extend(np.arange(b, e, s))
    times.append(end_time)
    return times


class endorse_2Dtest():

    def __init__(self, config):

        # TODO: set work dir
        self.work_dir = config["work_dir"]
        self.clean = config["clean_sample_dir"]
        self._config = config
        self.sample_dir = ""
        self.sample_counter = -1

    def set_parameters(self, data_par):
        param_list = self._config["surrDAMH_parameters"]["parameters"]
        assert(len(data_par) == len(param_list))

        for idx, param in enumerate(param_list):
            pname = param["name"]
            assert(pname in self._config["hm_params"])
            self._config["hm_params"][pname] = data_par[idx]

    def get_observations(self):
        try:
            print("get observations from flow_wrapper")
            res = self.calculate(self._config)
            return res
        except ValueError:
            print("flow_wrapper failed for unknown reason.")
            return -1000, []

    def calculate(self, config_dict):
        """
        The program changes to <work_dir> directory.
        does all the data preparation, passing
        running simulation
        extracting results
        """

        # create sample dir
        self.sample_counter = self.sample_counter + 1
        self.sample_dir = os.path.join(config_dict["work_dir"],
                                       "solver_" + str(config_dict["solver_id"]).zfill(2) +
                                       "_sample_" + str(self.sample_counter).zfill(3))
        os.makedirs(self.sample_dir, mode=0o775, exist_ok=True)
        os.chdir(self.sample_dir)

        print("=========================== RUNNING CALCULATION " +
              "solver {} ".format(config_dict["solver_id"]).zfill(2) +
              "sample {} ===========================".format(self.sample_counter).zfill(3),
              flush=True)
        print(self.sample_dir)

        # collect only
        if config_dict["collect_only"]:
            return 2, self.collect_results(config_dict)

        print("Creating mesh...")
        # comp_mesh = self.prepare_mesh(config_dict, cut_tunnel=False)
        comp_mesh = self.prepare_mesh(config_dict, cut_tunnel=True)

        mesh_bn = os.path.basename(comp_mesh)
        config_dict["hm_params"]["mesh"] = mesh_bn

        # endorse_2Dtest.read_physical_names(config_dict, comp_mesh)
        print("Creating mesh...finished")

        if config_dict["mesh_only"]:
            return -10, []  # tag, value_list

        # endorse_2Dtest.prepare_hm_input(config_dict)
        print("Running Flow123d - HM...")
        hm_succeed = self.call_flow(config_dict, 'hm_params', result_files=["flow_observe.yaml"])
        if not hm_succeed:
            # raise Exception("HM model failed.")
            # "Flow123d failed (wrong input or solver diverged)"
            print("Flow123d failed.")
            return -1, []  # tag, value_list
        print("Running Flow123d - HM...finished")

        if self._config["make_plots"]:
            try:
                self.observe_time_plot(config_dict)
            except:
                print("Making plot of sample results failed:")
                traceback.print_exc()
                return -2, []

        print("Finished computation")

        # collected_values = self.collect_results(config_dict)
        # print("Sample results collected.")
        # return 1, collected_values  # tag, value_list

        try:
            collected_values = self.collect_results(config_dict)
            print("Sample results collected.")
            return 1, collected_values  # tag, value_list
        except:
            print("Collecting sample results failed:")
            traceback.print_exc()
            return -3, []

    # def check_data(self, data, minimum, maximum):
    #     n_times = len(endorse_2Dtest.result_format()[0].times)
    #     if len(data) != n_times:
    #         raise Exception("Data not corresponding with time axis.")
    #
    #     if np.isnan(np.sum(data)):
    #         raise Exception("NaN present in extracted data.")
    #
    #     min = np.amin(data)
    #     if min < minimum:
    #         raise Exception("Data out of given range [min].")
    #     max = np.amax(data)
    #     if max > maximum:
    #         raise Exception("Data out of given range [max].")

    def collect_results(self, config_dict):
        output_dir = config_dict["hm_params"]["output_dir"]
        points2collect = config_dict["surrDAMH_parameters"]["observe_points"]

        # the times defined in input
        times = np.array(generate_time_axis(config_dict))
        with open(os.path.join(output_dir, "flow_observe.yaml"), "r") as f:
            loaded_yaml = yaml.load(f, yaml.CSafeLoader)
            points = loaded_yaml['points']
            point_names = [p["name"] for p in points]

            points2collect_indices = []
            for p2c in points2collect:
                tmp = [i for i, pn in enumerate(point_names) if pn == p2c]
                assert len(tmp) == 1
                points2collect_indices.append(tmp[0])

            print("Collecting results for observe points: ", points2collect)
            data = loaded_yaml['data']
            data_values = np.array([d["pressure_p0"] for d in data])
            values = data_values[:, points2collect_indices]
            obs_times = np.array([d["time"] for d in data]).transpose()

            # check that observe data are computed at all times of defined time axis
            all_times_computed = np.alltrue(np.isin(times, obs_times))
            if not all_times_computed:
                raise Exception("Observe data not computed at all times as defined by input!")
            # skip the times not specified in input
            t_indices = np.isin(obs_times, times).nonzero()
            values = values[t_indices].transpose()

        if config_dict["clean_sample_dir"]:
            shutil.rmtree(self.sample_dir)

        # flatten to format: [Point0_all_all_times, Point1_all_all_times, Point2_all_all_times, ...]
        res = values.flatten()
        return res

    def call_flow(self, config_dict, param_key, result_files):
        """
        Redirect sstdout and sterr, return true on succesfull run.
        :param result_files: Files to be computed - skip computation if already exist.
        :param param_key: config dict parameters key
        :param config_dict:
        :return:
        """

        status = False
        params = config_dict[param_key]
        fname = params["in_file"]
        # arguments = config_dict["_aux_flow_path"].copy()
        arguments = ['env', '-i']
        arguments.extend(config_dict["_aux_flow_path"].copy())
        output_dir = "output_" + fname
        config_dict[param_key]["output_dir"] = output_dir

        if all([os.path.isfile(os.path.join(output_dir, f)) for f in result_files]):
            status = True
        else:
            aux_functions.substitute_placeholders(
                os.path.join(config_dict["common_files_dir"],
                             fname + '_tmpl.yaml'),
                fname + '.yaml',
                params)

            arguments.extend(['--output_dir', output_dir, fname + ".yaml"])
            print("Running: ", " ".join(arguments))
            with open(fname + "_stdout", "w") as stdout:
                with open(fname + "_stderr", "w") as stderr:
                    completed = subprocess.run(arguments, stdout=stdout, stderr=stderr)
            print("Exit status: ", completed.returncode)
            status = completed.returncode == 0

        if status:
            log_file = os.path.join(self.sample_dir, output_dir, "flow123.0.log")
            conv_check = aux_functions.check_conv_reasons(log_file)
            print("converged: ", conv_check)
            status = conv_check >= 0

        return status


    def prepare_mesh(self, config_dict, cut_tunnel):
        mesh_name = config_dict["geometry"]["mesh_name"]
        if cut_tunnel:
            mesh_name = mesh_name + "_cut"
        mesh_file = mesh_name + ".msh"
        mesh_healed = mesh_name + "_healed.msh"

        # suppose that the mesh was created/copied during preprocess
        assert os.path.isfile(os.path.join(config_dict["common_files_dir"], mesh_healed))
        shutil.copyfile(os.path.join(config_dict["common_files_dir"], mesh_healed), mesh_healed)
        return mesh_healed


    def observe_time_plot(self, config_dict):

        output_dir = config_dict["hm_params"]["output_dir"]

        with open(os.path.join(output_dir, "flow_observe.yaml"), "r") as f:
            loaded_yaml = yaml.load(f, yaml.CSafeLoader)
            points = loaded_yaml['points']
            point_names = [p["name"] for p in points]
            data = loaded_yaml['data']
            values = np.array([d["pressure_p0"] for d in data]).transpose()
            times = np.array([d["time"] for d in data]).transpose()

            fig, ax1 = plt.subplots()
            temp_color = ['red', 'green', 'violet', 'blue']
            ax1.set_xlabel('time [d]')
            ax1.set_ylabel('pressure [m]')
            for i in range(0, len(point_names)):
                ax1.plot(times, values[i, 0:], color=temp_color[i], label=point_names[i])

            ax1.tick_params(axis='y')
            ax1.legend()

            fig.tight_layout()  # otherwise the right y-label is slightly clipped
            # plt.show()
            plt.savefig("observe_pressure.pdf")