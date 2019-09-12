import os
import sys
import numpy as np
import threading
import subprocess
import yaml
import attr
from typing import Any, List, Dict, Tuple

src_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(src_path, '../MLMC/src'))
sys.path.append(os.path.join(src_path, '../dfn/src'))
from mlmc import flow_mc
from geomop import polygons
import mlmc.random.fracture as fracture


def in_file(base):
    return "flow_{}.yaml".format(base)

def mesh_file(base):
    return "mesh_{}.msh".format(base)

def fields_file(base):
    return "fields_{}.msh".format(base)



class FlowThread(threading.Thread):

    def __init__(self, basename, config_dict):
        self.base = basename
        self.flow_args = config_dict["flow_executable"].copy()
        super().__init__()

    def run(self):
        in_f = in_file(self.base)
        out_dir = self.base

        flow_in = "flow_{}.yaml".format(run_basename)
        flow_mc.substitute_placeholders("flow_templ.yaml", in_f,
                                        mesh_file = mesh_file(self.base),
                                        fields_file=fields_file(self.base))
        self.flow_args.extend(['--output_dir', out_dir, in_f])
        with open(self.base + "_stdout", "w") as stdout:
            with open(self.base + "_stderr", "w") as stderr:
                completed = subprocess.run(self.flow_args, stdout=stdout, stderr=stderr)
            print("Exit status: ", completed.returncode)
            status = completed.returncode == 0
        conv_check = self.check_conv_reasons(os.path.join(out_dir, "flow123.0.log"))
        print("converged: ", conv_check)
        return status  # and conv_check

    def check_conv_reasons(self, log_fname):
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


@attr.s(auto_attribs=True)
class Region:
    name:str = ""
    dim:int = -1
    boundary:bool = False
    mesh_step:float = 0.0


@attr.s(auto_attribs=True)
class FlowProblem:
    basename: str
    fr_range: Tuple[float, float]
    fractures: List[Any]
    config_dict: Dict[str, Any]


    def make_fracture_network(self):
        from geomop.plot_polygons import plot_decomp_segments
        self.mesh_step = self.fr_range[0]

        # Init regions
        self.regions={}
        self.regions['none'] = none_reg = Region('none', 0)
        self.regions['bulk'] = Region('bulk', dim=2, mesh_step=self.mesh_step)


        pd = polygons.PolygonDecomposition([none_reg, none_reg, none_reg], self.mesh_step)

        # make outer polygon
        geom = self.config_dict["geometry"]
        lx, ly = geom["domain_box"]
        outer_polygon = [[-lx / 2, -ly / 2], [+lx / 2, -ly / 2], [+lx / 2, +ly / 2], [-lx / 2, +ly / 2]]
        last_pt = outer_polygon[-1]
        side_segments = {}
        for i_side, pt in enumerate(outer_polygon):
            reg = Region(".side_{}".format(i_side), dim=1, mesh_step=self.mesh_step)
            self.regions[reg.name] = reg
            sub_segments = pd.add_line(last_pt, pt, attr=reg)
            last_pt = pt
            assert type(sub_segments) == list and len(sub_segments) == 1
            seg = sub_segments[0]
            side_segments[seg.id] = i_side
        assert len(pd.polygons) == 2
        pd.polygons[1].attr = self.regions['bulk']

        # extract fracture lines larger then the mesh step
        # self.fracture_lines = self.fractures.get_lines(self.fr_range)
        # outer_wire = pd.outer_polygon.outer_wire.childs
        # assert len(outer_wire) == 1
        # outer_wire = next(iter(outer_wire))
        # for i_fr, (p0, p1) in self.fracture_lines.items():
        #     reg = Region("fr_{}".format(i_fr), dim=1, mesh_step=self.mesh_step)
        #     self.regions[reg.name] = reg
        #     print(i_fr, "fr size:", np.linalg.norm(p1 - p0))
        #     try:
        #         segments = pd.add_line(p0, p1, attr=reg)
        #         w = pd.decomp.wires.get(771, None)
        #         if w is not None and (w.parent == w or w in w.childs):
        #             assert False
        #     except Exception as e:
        #         # new_points = [pt for seg in segments for pt in seg.vtxs]
        #         # plot_decomp_segments(pd, [p0, p1])
        #         # raise e
        #         # print("False")
        #         pass
        #     # pd.decomp.check_consistency()
        #
        #     # remove segments out of the outer polygon
        #     if type(segments) == list:
        #         for seg in segments:
        #             if seg.wire[0] == seg.wire[1] and seg.wire[0] == outer_wire:
        #                 points = seg.vtxs
        #                 pd.delete_segment(seg)
        #                 for pt in points:
        #                     if pt.is_free():
        #                         pd.remove_free_point(pt.id)
        #
        # # plot_decomp_segments(pd)
        self.decomp = pd

    def make_mesh(self):
        import geometry_2d as geom

        gmsh_executable = self.config_dict["gmsh_executable"]
        g2d = geom.Geometry2d("fine", self.regions)
        g2d.add_compoud(self.decomp)
        g2d.make_brep_geometry()
        step_range = (self.mesh_step * 0.9, self.mesh_step *1.1)
        g2d.call_gmsh(gmsh_executable, step_range)
        self.mesh = g2d.modify_mesh()

    def make_fields(self):
        pass

    def run(self):
        self.thread = FlowThread(self.basename, self.config_dict)
        return self.thread

    def elementwise_mesh(self, coarse_mesh):
        pass

class BothSample:


    def __init__(self, sample_config):
        """
        Attributes:
        seed,
        h_fine_step,
        h_coarse_step,
        do_coarse
        config_path
        :param sample_config:
        """
        self.__dict__.update(sample_config)
        np.random.seed(self.seed)
        with open(self.config_path, "r") as f:
            self.config_dict = yaml.load(f, Loader=yaml.FullLoader)

    def generate_fractures(self):
        geom = self.config_dict["geometry"]
        lx, ly = geom["fractures_box"]
        fr_size_range = geom["pow_law_size_range"]
        pow_law_exp_3d = geom["pow_law_size_exp"]
        pow_law_sample_range = geom["pow_law_sample_range"]
        n_frac_limit = geom["n_frac_limit"]
        p_32 = geom["p_32"]

        # generate fracture set
        fracture_box = [lx, ly, 0]
        area = lx * ly

        pop = fracture.Population(area, fracture.LineShape)
        pop.add_family("all",
                       fracture.FisherOrientation(0, 90, np.inf),
                       fracture.VonMisesOrientation(0, 0),
                       fracture.PowerLawSize.from_mean_area(pow_law_exp_3d - 1, fr_size_range, p_32, pow_law_exp_3d))

        if pow_law_sample_range:
            pop.set_sample_range(pow_law_sample_range)
        elif n_frac_limit:
            pop.set_sample_range([None, max(lx, ly)], sample_size=n_frac_limit)

        print("total mean size: ", pop.mean_size())
        print("size range:", pop.families[0].size.sample_range)
        pos_gen = fracture.UniformBoxPosition(fracture_box)
        fractures = pop.sample(pos_distr=pos_gen, keep_nonempty=True)

        fr_set = fracture.Fractures(fractures, fr_size_range[0] / 2)
        return fr_set


    def calculate(self):
        fractures = self.generate_fractures()
        # fine problem
        fine_flow = FlowProblem("fine", (self.h_fine_step, np.inf), fractures, self.config_dict)
        fine_flow.make_fracture_network()
        fine_flow.make_mesh()
        fine_flow.make_fields()
        fine_flow.run()

        # coarse problem
        if self.do_coarse:
            coarse_flow = FlowProblem("coarse", (self.h_coarse_step, np.inf), fractures, self.config_dict)
            coarse_flow.make_fracture_network()
            coarse_flow.make_mesh()

            coarse_ref = FlowProblem("coarse_ref", (0, self.h_coarse_step), fractures, self.config_dict)
            coarse_ref.make_fracture_network()
            coarse_ref.elementwise_mesh(coarse_flow.mesh)
            coarse_ref.make_fields()
            coarse_ref.run().join()
            cond = coarse_ref.element_conductivities()

            coarse_flow.make_fields(cond)
            coarse_flow.run().join()
        fine_flow.thread.join()







    # def mean_tensor(self, balance_dict, normals, regions):
    #     for n, reg in zip(normals, regions):
    #         balance_dict





if __name__ == "__main__":
    sample_config = sys.argv[1]
    with open(sample_config, "r") as f:
        config_dict = yaml.load(f, Loader=yaml.FullLoader)
    bs = BothSample(config_dict)
    bs.calculate()

