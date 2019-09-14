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

    def __init__(self, basename, outer_regions, config_dict):
        self.base = basename
        self.outer_regions_list = outer_regions
        self.flow_args = config_dict["flow_executable"].copy()
        n_steps = 30
        t = np.linspace(0.0, np.pi, n_steps)
        self.p_loads = np.array([np.cos(t), np.sin(t)]).T
        super().__init__()

    def run(self):
        in_f = in_file(self.base)
        out_dir = self.base
        n_loads = len(self.p_loads)
        flow_in = "flow_{}.yaml".format(self.base)
        params = dict(
            mesh_file=mesh_file(self.base),
            fields_file=fields_file(self.base),
            outer_regions=str(self.outer_regions_list),
            n_steps=len(self.p_loads)
            )
        flow_mc.substitute_placeholders("flow_templ.yaml", in_f, params)
        self.flow_args.extend(['--output_dir', out_dir, in_f])

        if os.path.exists(os.path.join(out_dir, "water_balance.yaml")):
            return True
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

    def is_active(self, dim):
        active = dim >= 0
        if active:
            assert dim == self.dim, "Can not create shape of dim: {} in region '{}' of dim: {}.".format(dim, self.name, self.dim)
        return active




@attr.s(auto_attribs=True)
class FlowProblem:
    basename: str
    fr_range: Tuple[float, float]
    fractures: List[Any]
    config_dict: Dict[str, Any]

    regions: List[Region] = []
    side_regions: List[Region] = []
    reg_to_fr: Dict[int, int] = {}

    def add_region(self, name, dim, mesh_step=0.0, boundary=False):
        reg = Region(name, dim, boundary, mesh_step)
        reg.id = len(self.regions)
        self.regions.append(reg)
        return reg


    def make_fracture_network(self):
        from geomop.plot_polygons import plot_decomp_segments
        self.mesh_step = self.fr_range[0]

        # Init regions
        none_reg = self.add_region('none', dim=-1)
        bulk_reg = self.add_region('bulk', dim=2, mesh_step=self.mesh_step)

        pd = polygons.PolygonDecomposition(self.mesh_step)
        # make outer polygon
        geom = self.config_dict["geometry"]
        lx, ly = geom["domain_box"]
        outer_polygon = [[-lx / 2, -ly / 2], [+lx / 2, -ly / 2], [+lx / 2, +ly / 2], [-lx / 2, +ly / 2]]
        last_pt = outer_polygon[-1]
        side_segments = {}
        for i_side, pt in enumerate(outer_polygon):
            reg = self.add_region(".side_{}".format(i_side), dim=1, mesh_step=self.mesh_step, boundary=True)
            reg.sub_reg = self.add_region(".side_fr_{}".format(i_side), dim=0, mesh_step=self.mesh_step, boundary=True)
            self.side_regions.append(reg)
            diff = np.array(pt) - np.array(last_pt)
            normal = np.array([diff[1], -diff[0]])
            reg.normal =  normal / np.linalg.norm(normal)
            sub_segments = pd.add_line(last_pt, pt)
            for seg in sub_segments:
                seg.attr = reg
            last_pt = pt
            assert type(sub_segments) == list and len(sub_segments) == 1
            seg = sub_segments[0]
            side_segments[seg.id] = i_side
        assert len(pd.polygons) == 2
        pd.polygons[1].attr = bulk_reg

        # extract fracture lines larger then the mesh step
        self.fracture_lines = self.fractures.get_lines(self.fr_range)
        self.fracture_lines = {0: [np.array([0, -50]), np.array([0, 50])]}
        outer_wire = pd.outer_polygon.outer_wire.childs
        assert len(outer_wire) == 1
        outer_wire = next(iter(outer_wire))
        for i_fr, (p0, p1) in self.fracture_lines.items():
            reg = self.add_region("fr_{}".format(i_fr), dim=1, mesh_step=self.mesh_step)
            self.reg_to_fr[reg.id] = i_fr
            print(i_fr, "fr size:", np.linalg.norm(p1 - p0))
            try:
                sub_segments = pd.add_line(p0, p1)
            except Exception as e:
                # new_points = [pt for seg in segments for pt in seg.vtxs]
                # plot_decomp_segments(pd, [p0, p1])
                print(e)
                pass
            # pd.decomp.check_consistency()

            # remove segments out of the outer polygon
            if type(sub_segments) == list:
                for seg in sub_segments:
                    if seg.attr is None:
                        seg.attr = reg
                    if seg.wire[0] == seg.wire[1] and seg.wire[0] == outer_wire:
                        points = seg.vtxs
                        pd.delete_segment(seg)
                        for pt in points:
                            if pt.is_free():
                                pd.remove_free_point(pt.id)

        # assign boundary region to outer polygion points
        for seg, side in outer_wire.segments():
            side_reg = seg.attr
            assert hasattr(side_reg, "sub_reg")
            seg.vtxs[side].attr = side_reg.sub_reg
        # none region to remaining
        for shape_list in pd.decomp.shapes:
            for shape in shape_list.values():
                if shape.attr is None:
                    shape.attr = none_reg

        # # plot_decomp_segments(pd)
        self.decomp = pd

    def make_mesh(self):
        import geometry_2d as geom

        gmsh_executable = self.config_dict["gmsh_executable"]
        g2d = geom.Geometry2d("mesh_" + self.basename, self.regions)
        g2d.add_compoud(self.decomp)
        g2d.make_brep_geometry()
        step_range = (self.mesh_step * 0.9, self.mesh_step *1.1)
        g2d.call_gmsh(gmsh_executable, step_range)
        self.mesh = g2d.modify_mesh()

    def make_fields(self):
        bulk_elements = [eid
                         for eid, (t, tags, nodes) in self.mesh.elements.items()
                         if not self.regions[tags[0]-10000].boundary]
        I_tn = np.eye(3, dtype=float)
        const_conductivity = float(self.config_dict['bulk_conductivity']) * I_tn
        aperture_per_size = float(self.config_dict['aperture_per_size'])
        n_elem = len(bulk_elements)
        elem_ids = []
        cond_tn_field = np.empty((n_elem, 9))
        cs_field = np.empty((n_elem, 1))
        for el_id in bulk_elements:
            ele = self.mesh.elements[el_id]
            el_type, tags, node_ids = ele
            n_nodes = len(node_ids)
            if n_nodes == 2:
                reg_id = tags[0] - 10000
                # line, compute conductivity from fracture size using cubic law
                i_fr = self.reg_to_fr[reg_id]
                fr_size = self.fractures.fractures[i_fr].rx
                cs = fr_size * aperture_per_size
                cond = cs ** 2 / 12
                cond_tn = cond * I_tn
            else:
                cs = 1.0
                cond_tn = const_conductivity

            i = len(elem_ids)
            elem_ids.append(el_id)
            cond_tn_field[i] = cond_tn.flatten()
            cs_field[i] = cs


        fname = fields_file(self.basename)
        with open(fname, "w") as fout:
            self.mesh.write_ascii(fout)
            self.mesh.write_element_data(fout, elem_ids, 'conductivity_tensor', cond_tn_field)
            self.mesh.write_element_data(fout, elem_ids, 'cross_section', cs_field)


    def elementwise_mesh(self, coarse_mesh):
        import geometry_2d as geom

        self.reg_to_coarse_el = {}  # bulk and fracture region id to coarse element id
        g2d = geom.Geometry2d("mesh_" + self.basename, self.regions)
        for eid, (tele, tags, nodes) in coarse_mesh.elements.items():
            # create regions
            # outer polygon
            # add fractures
            # add compound, mark regions

        g2d.add_compoud(self.decomp)
        g2d.make_brep_geometry()

        step_range = (self.mesh_step * 0.9, self.mesh_step *1.1)
        gmsh_executable = self.config_dict["gmsh_executable"]
        g2d.call_gmsh(gmsh_executable, step_range)
        self.mesh = g2d.modify_mesh()


    def run(self):
        outer_reg_names = []
        for reg in self.side_regions:
            outer_reg_names.append(reg.name)
            outer_reg_names.append(reg.sub_reg.name)
        self.thread = FlowThread(self.basename, outer_reg_names, self.config_dict)

        self.thread.start()
        return self.thread

    def effective_tensor(self, side_regions):
        """
        :param mesh: GmshIO mesh object.
        :param side_regions: List of side regions with the "normal" attribute.
        :return:
        """
        loads = self.thread.p_loads
        with open(os.path.join(self.basename, "water_balance.yaml")) as f:
            balance = yaml.load(f, Loader=yaml.FullLoader)['data']
        flux_response = np.zeros_like(loads)
        reg_map = {}
        for reg in side_regions:
            reg_map[reg.name] = reg
            reg_map[reg.sub_reg.name] = reg
        for entry in balance:
            reg = reg_map.get(entry['region'], None)
            bc_influx = entry['data'][0]
            if reg is not None:
                flux_response[entry['time']] += reg.normal * bc_influx
        flux_response /= len(side_regions)
        #flux_response *= np.array([100, 1])[None, :]


        # least square fit for the symmetric conductivity tensor
        rhs = flux_response.flatten()
        # columns for the tensor values: C00, C01, C11
        pressure_matrix = np.zeros((len(rhs), 3))
        for i_load, (p0, p1) in enumerate(loads):
            i0 = 2 * i_load
            i1 = i0 + 1
            pressure_matrix[i0] = [p0, p1, 0]
            pressure_matrix[i1] = [0, p0, p1]
        C = np.linalg.lstsq(pressure_matrix, rhs, rcond=None)[0]
        cond_tn = np.array([[C[0], C[1]], [C[1], C[2]]])
        self.plot_effective_tensor(loads, flux_response, cond_tn)
        print(cond_tn)
        return cond_tn

    def plot_effective_tensor(self, loads, fluxes, cond_tn):
        """
        Plot response fluxes for pressure gradients with angles [0, pi) measured from the X-azis counterclockwise.
        :param fluxes: np.array of 2d fluex vectors.
        :return:
        """

        import matplotlib.pyplot as plt

        fig = plt.figure()
        # setting the axis limits in [left, bottom, width, height]
        #rect = [0.1, 0.1, 0.8, 0.8]
        ax = fig.add_subplot()

        #ax_polar = fig.add_axes(rect, polar=True, frameon=False)
        X, Y = cond_tn @ loads.T
        ax.scatter(X, Y, c='green')
        ax.scatter(-X, -Y, c='green')
        ax.scatter(fluxes[:, 0], fluxes[:, 1], c='red')
        ax.scatter(-fluxes[:, 0], -fluxes[:, 1], c='red')
        ax.set_xlim(-max(X), max(X))
        ax.set_ylim(-max(Y), max(Y))

        #ax_polar.grid(True)
        fig.savefig("conductivity_tensor.pdf")
        plt.show()





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
            #coarse_ref.make_fracture_network()
            coarse_ref.elementwise_mesh(coarse_flow.mesh)
            coarse_ref.make_fields()
            coarse_ref.run().join()
            cond = coarse_ref.element_conductivities()

            coarse_flow.make_fields(cond)
            coarse_flow.run().join()
        fine_flow.thread.join()
        fine_flow.effective_tensor(fine_flow.side_regions)







    # def mean_tensor(self, balance_dict, normals, regions):
    #     for n, reg in zip(normals, regions):
    #         balance_dict





if __name__ == "__main__":
    sample_config = sys.argv[1]
    with open(sample_config, "r") as f:
        config_dict = yaml.load(f, Loader=yaml.FullLoader)
    bs = BothSample(config_dict)
    bs.calculate()

