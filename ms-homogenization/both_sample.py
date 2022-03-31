from typing import *
from abc import *

import os
import sys
import numpy as np
import threading
import subprocess
import yaml
import attr
import collections
import traceback
import time
import pandas
import scipy.spatial as sc_spatial
import scipy.interpolate as sc_interpolate
import atexit
from scipy import stats
from scipy.spatial import distance

src_path = os.path.dirname(os.path.abspath(__file__))

from bgem.gmsh import gmsh_io
from bgem.polygons import polygons
import fracture
import matplotlib.pyplot as plt
import copy
import skgstat


def in_file(base):
    return "flow_{}.yaml".format(base)


def mesh_file(base):
    return "mesh_{}.msh".format(base)


def fields_file(base):
    return "fields_{}.msh".format(base)


def substitute_placeholders(file_in, file_out, params):
    """
    Substitute for placeholders of format '<name>' from the dict 'params'.
    :param file_in: Template file.
    :param file_out: Values substituted.
    :param params: { 'name': value, ...}
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


class FlowThread(threading.Thread):

    def __init__(self, basename, outer_regions, config_dict):
        self.base = basename
        self.outer_regions_list = outer_regions
        self.flow_args = config_dict["flow_executable"].copy()
        n_steps = config_dict["n_pressure_loads"]
        t = np.pi * np.arange(0, n_steps) / n_steps
        self.p_loads = np.array([np.cos(t), np.sin(t)]).T
        super().__init__()

    def run(self):
        in_f = in_file(self.base)
        out_dir = self.base
        # n_loads = len(self.p_loads)
        # flow_in = "flow_{}.yaml".format(self.base)
        params = dict(
            mesh_file=mesh_file(self.base),
            fields_file=fields_file(self.base),
            outer_regions=str(self.outer_regions_list),
            n_steps=len(self.p_loads)
            )
        substitute_placeholders("flow_templ.yaml", in_f, params)
        self.flow_args.extend(['--output_dir', out_dir, in_f])

        if os.path.exists(os.path.join(out_dir, "flow_fields.msh")):
            return True
        with open(self.base + "_stdout", "w") as stdout:
            with open(self.base + "_stderr", "w") as stderr:
                print("flow args ", self.flow_args)
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
        active = self.dim >= 0
        if active:
            assert dim == self.dim, "Can not create shape of dim: {} in region '{}' of dim: {}.".format(dim, self.name, self.dim)
        return active


def gmsh_mesh_bulk_elements(mesh):
    """
    Generator of IDs of bulk elements.
    :param mesh:
    :return:
    """
    is_bc_reg_id={}
    for name, reg_id_dim in mesh.physical.items():
        is_bc_reg_id[reg_id_dim] = (name[0] == '.')
    for eid, ele in mesh.elements.items():
        (t, tags, nodes) = ele
        dim = len(nodes) - 1
        if not is_bc_reg_id[(tags[0], dim)]:
            yield eid, ele


class BulkBase(ABC):
    @abstractmethod
    def element_data(self, mesh, eid):
        """
        :return:
        """
        pass


@attr.s(auto_attribs=True)
class BulkFields(BulkBase):
    mean_log_conductivity: Tuple[float, float]
    cov_log_conductivity: Optional[List[List[float]]]
    angle_mean: float = attr.ib(converter=float)
    angle_concentration: float = attr.ib(converter=float)

    def element_data(self, mesh, eid):
        # Unrotated tensor (eigenvalues)
        if self.cov_log_conductivity is None:
            log_eigenvals = self.mean_log_conductivity
        else:
            log_eigenvals = np.random.multivariate_normal(
                mean=self.mean_log_conductivity,
                cov=self.cov_log_conductivity,
                )

        unrotated_tn = np.diag(np.power(10, log_eigenvals))

        # rotation angle
        if self.angle_concentration is None or self.angle_concentration == 0:
            angle = np.random.uniform(0, 2*np.pi)
        elif self.angle_concentration == np.inf:
            angle = self.angle_mean
        else:
            angle = np.random.vonmises(self.angle_mean, self.angle_concentration)
        c, s = np.cos(angle), np.sin(angle)
        rot_mat = np.array([[c, -s], [s, c]])
        cond_2d = rot_mat @ unrotated_tn @ rot_mat.T
        return 1.0, cond_2d


class BulkMicroScale(BulkBase):
    def __init__(self, microscale):
        self.microscale = microscale
        self.microscale_tensors = None

    def element_data(self, mesh, eid):
        if self.microscale_tensors is None:
            self.microscale_tensors = self.microscale.effective_tensor_from_bulk()
        return 1.0, self.microscale_tensors[eid]


class BulkFromFine(BulkBase):
    def __init__(self, fine_problem):
        points, values = fine_problem.bulk_field()
        self.mean_val = np.mean(values)
        tria = sc_spatial.Delaunay(points)
        #print("Values, shape:", values.shape)
        self.interp =  sc_interpolate.LinearNDInterpolator(tria, values.T, fill_value=0)
        self.interp_nearest = sc_interpolate.LinearNDInterpolator(points, values.T)

    def element_data(self, mesh, eid):
        el_type, tags, node_ids = mesh.elements[eid]
        center = np.mean([np.array(mesh.nodes[nid]) for nid in node_ids], axis=0)
        v = self.interp(center[0:2])
        if np.alltrue(v == 0):
            #v = self.interp_nearest(center[0:2])
            v = [[self.mean_val,0,self.mean_val]]
        #print("V, shape:", v.shape)
        v00, v01, v11 = v[0]
        cond = np.array([[v00, v01], [v01, v11]])
        #print("cond, shape: ", cond.shape)
        e0, e1 = np.linalg.eigvalsh(cond)
        if e0 < 1e-20 or e1 < 1e-20:
            print(e0, e1, v, center)
            assert False
        return 1.0, cond


class BulkChoose(BulkBase):
    def __init__(self, finer_level_path):
        self.cond_tn = np.array(pandas.read_csv(finer_level_path, sep=' '))

    def element_data(self, mesh, eid):
        idx = np.random.randint(len(self.cond_tn))
        return 1.0, self.cond_tn[idx].reshape(2,2)


@attr.s(auto_attribs=True)
class FractureModel:
    fractures: fracture.Fractures
    region_to_fracture: Dict[int, int]
    aperture_per_size: float = attr.ib(converter=float)
    water_viscosity: float = attr.ib(converter=float)
    gravity_accel: float = attr.ib(converter=float)
    water_density: float = attr.ib(converter=float)

    def element_data(self, mesh, eid):
        el_type, tags, node_ids = mesh.elements[eid]
        reg_id = tags[0] - 10000
        # line, compute conductivity from fracture size using cubic law
        # Isotropic conductivity in fractures. (Simplification.)
        i_fr = self.region_to_fracture[reg_id]
        fr_size = self.fractures.fractures[i_fr].rx
        cs = fr_size * self.aperture_per_size
        cond = cs ** 2 / 12 * self.water_density * self.gravity_accel / self.water_viscosity
        # print(f"fr: {fr_size} {cs} {cond}")
        cond_tn = cond * np.eye(2, 2)
        return cs, cond_tn


def tensor_3d_flatten(tn_2d):
    tn3d = np.eye(3)
    tn3d[0:2, 0:2] = tn_2d
    # tn3d[0:2, 0:2] += tn_2d # ???
    return tn3d.ravel()


def write_fields(mesh, basename, bulk_model, fracture_model):
    elem_ids = []
    cond_tn_field = []
    cs_field = []
    print("mesh ", mesh)
    for el_id, ele in gmsh_mesh_bulk_elements(mesh):
        # print("ele id ", el_id)
        # print("ele ", ele)
        elem_ids.append(el_id)
        el_type, tags, node_ids = ele
        n_nodes = len(node_ids)

        if n_nodes == 2:
            cs, cond_tn = fracture_model.element_data(mesh, el_id)
        else:
            # n_nodes == 3
            cs, cond_tn = bulk_model.element_data(mesh, el_id)

        cs_field.append(np.array(cs))
        cond_tn_field.append(tensor_3d_flatten(cond_tn))

    fname = fields_file(basename)
    with open(fname, "w") as fout:
        mesh.write_ascii(fout)

        # print("fout ", fout)
        # print("elem_ids ", elem_ids)
        # print("cond tn field ", cond_tn_field)
        mesh.write_element_data(fout, elem_ids, 'conductivity_tensor', np.array(cond_tn_field))
        mesh.write_element_data(fout, elem_ids, 'cross_section', np.array(cs_field).reshape(-1, 1))
    return elem_ids, cs_field, cond_tn_field


@attr.s(auto_attribs=True)
class FlowProblem:
    i_level: int
    # MLMC Level index (to retrieve model parameters from main config)
    basename: str
    # Basename for files of this flow problem.
    fr_range: Tuple[float, float]
    # Fracture range to extract from the full list of the generated fractures.
    fractures: List[Any]
    # The Fractures object with generated fractures.
    bulk_model: BulkBase
    # The bulk model (specific for fine, coarse, etc.)
    config_dict: Dict[str, Any]
    # global config dictionary.

    regions: List[Region] = attr.ib(factory=list)
    # List of regions used in the geometry and mesh preparation
    side_regions: List[Region] = attr.ib(factory=list)
    # separate regions for sides of the outer wire, with fracture subregions and normals
    # used for creating boundary region set and for boundary averaging of the conductivity tensor
    reg_to_fr: Dict[int, int] = attr.ib(factory=dict)
    # Maps region id to original fracture id
    reg_to_group: List[Tuple[int, int]] = attr.ib(factory=dict)
    # Groups of regions for which the effective conductivity tensor will be computed separately.
    # One group is specified by the tuple of bulk and fracture region ID.
    group_positions: Dict[int, np.array] = attr.ib(factory=dict)
    # Centers of macro elements.
    skip_decomposition: bool = False

    # created later
    mesh: gmsh_io.GmshIO = None

    # safe conductivities produced by `make_fields`
    _elem_ids: Any = None
    _cond_tn_field: Any = None

    @classmethod
    def make_fine(cls, i_level, fr_range, fractures, finer_level_path, config_dict):
        level_dict = config_dict['levels'][i_level]
        bulk_conductivity = level_dict['bulk_conductivity']
        print("bulk conductivity ", bulk_conductivity)

        if bulk_conductivity.get('choose_from_finer_level', False):
            bulk_model = BulkChoose(finer_level_path)
        else:
            bulk_model = BulkFields(**bulk_conductivity)
        return FlowProblem(i_level, "fine",
                           fr_range, fractures, bulk_model, config_dict)

    @classmethod
    def make_coarse(cls, i_level, fr_range, fractures, micro_scale_problem, config_dict):
        bulk_model = BulkMicroScale(micro_scale_problem)
        return FlowProblem(i_level, "coarse",
                           fr_range, fractures, bulk_model, config_dict)

    @classmethod
    def make_microscale(cls, i_level, fr_range, fractures, fine_flow, config_dict):
        # use bulk fields from the fine level
        bulk_model = BulkFromFine(fine_flow)

        return FlowProblem(i_level, "coarse_ref",
                           fr_range, fractures, bulk_model, config_dict)

    @property
    def pressure_loads(self):
        return self.thread.p_loads

    def add_region(self, name, dim, mesh_step=0.0, boundary=False):
        reg = Region(name, dim, boundary, mesh_step)
        reg.id = len(self.regions)
        self.regions.append(reg)
        return reg

    def init_decomposition(self, outer_polygon, bulk_reg, tol):
        """
        Create polygon decomposition and add the outer polygon.
        :param outer_polygon: [np.array[2], ..] Vertices of the outer polygon.
        :return: (PolygonDecomposition, side_regions)
        """
        pd = polygons.PolygonDecomposition(tol)
        last_pt = outer_polygon[-1]
        print("last pt ", last_pt)

        side_regions = []
        for i_side, pt in enumerate(outer_polygon):
            reg = self.add_region(".side_{}".format(i_side), dim=1, mesh_step=self.mesh_step, boundary=True)
            print("reg ", reg)
            print("i side ", i_side)

            reg.sub_reg = self.add_region(".side_fr_{}".format(i_side), dim=0, mesh_step=self.mesh_step, boundary=True)
            diff = np.array(pt) - np.array(last_pt)
            normal = np.array([diff[1], -diff[0]])
            reg.normal = normal / np.linalg.norm(normal)
            side_regions.append(reg)

            # print("last pt ", last_pt)
            # print("pt ", pt)

            sub_segments = pd.add_line(last_pt, pt, deformability=0)

            #print("sub segments ", sub_segments)

            if not (type(sub_segments) == list and len(sub_segments) == 1):
                from bgem.polygons.plot_polygons import plot_decomp_segments
                plot_decomp_segments(pd)

            for seg in sub_segments:
                seg.attr = reg
            last_pt = pt
            #assert type(sub_segments) == list and len(sub_segments) == 1, sub_segments
        assert len(pd.polygons) == 2
        pd.polygons[1].attr = bulk_reg

        print("side regions ", side_regions)

        print("pd.points ", pd.points)

        return pd, side_regions

    def add_fractures(self, pd, fracture_lines, eid):
        outer_wire = pd.outer_polygon.outer_wire.childs
        #print("outer wire ", outer_wire)

        assert len(outer_wire) == 1
        outer_wire = next(iter(outer_wire))
        fracture_regions = []
        for i_fr, (p0, p1) in fracture_lines.items():
            reg = self.add_region("fr_{}".format(i_fr), dim=1, mesh_step=self.mesh_step)
            self.reg_to_fr[reg.id] = i_fr
            fracture_regions.append(reg)
            if self.skip_decomposition:
                continue
            #print("    ", i_fr, "fr size:", np.linalg.norm(p1 - p0))
            try:
                #pd.decomp.check_consistency()
                # if eid == 0 and i_fr == 5:
                #      print("stop")
                #      # plot_decomp_segments(pd, [p0, p1])
                # TODO: make short fractures more deformable
                sub_segments = pd.add_line(p0, p1, deformability=1)
            except Exception as e:
                # new_points = [pt for seg in segments for pt in seg.vtxs]
                print('Decomp Error, dir: {} base: {} eid: {}  i_fr: {}'.format(os.getcwd(), self.basename, eid, i_fr))
                traceback.print_exc()
                #plot_decomp_segments(pd, [p0, p1])
                #raise e
                #print(e)
                #pass
            # pd.decomp.check_consistency()

            # remove segments out of the outer polygon
            if type(sub_segments) == list:
                #print("sub segments ", sub_segments)
                #exit()
                for seg in sub_segments:
                    if seg.attr is None:
                        seg.attr = reg
                    # remove segments poking out of outer polygon
                    if seg.wire[0] == seg.wire[1] and (seg.wire[0].parent == pd.outer_polygon.outer_wire):
                        points = seg.vtxs
                        # print("pd segments ", pd.segments.keys())
                        # print("seg ", seg)
                        if seg.id in pd.segments.keys():
                            pd.delete_segment(seg)
                            for pt in points:
                                if pt.is_free():# and pt.id in pd.points.keys():
                                    pd._rm_point(pt)
                                # else:
                                #     print("not rm pt.id ", pt.id)

        #plot_decomp_segments(pd, [p0, p1])
        # assign boundary region to outer polygon points
        for seg, side in outer_wire.segments():
            side_reg = seg.attr
            if side_reg.boundary:
                assert hasattr(side_reg, "sub_reg")
                seg.vtxs[side].attr = side_reg.sub_reg
        # none region to remaining
        for shape_list in pd.decomp.shapes:
            for shape in shape_list.values():
                if shape.attr is None:
                    shape.attr = self.none_reg

        # # plot_decomp_segments(pd)
        return pd, fracture_regions

    def make_fracture_network(self):
        self.mesh_step = self.fr_range[0]

        print("self. fr range ", self.fr_range)

        # Init regions
        self.none_reg = self.add_region('none', dim=-1)
        bulk_reg = self.add_region('bulk_2d', dim=2, mesh_step=self.mesh_step)
        self.reg_to_group[bulk_reg.id] = 0

        # make outer polygon
        geom = self.config_dict["geometry"]
        lx, ly = geom["domain_box"]

        if "outer_polygon" in geom:
            self.outer_polygon = geom["outer_polygon"]
        else:
            self.outer_polygon = [[-lx / 2, -ly / 2], [+lx / 2, -ly / 2], [+lx / 2, +ly / 2], [-lx / 2, +ly / 2]]

        #self.outer_polygon = [[-lx / 2, 0], [0, 0], [0, +ly / 2], [-lx / 2, +ly / 2]]

        #self.outer_polygon = [[-10, -10], [10, -10], [10, 10], [-10, 10]]
        #self.outer_polygon = [[-500, 0], [0, 0], [0, 500], [-500, 500]]

        print("self outer polygon ", self.outer_polygon)

        pd, self.side_regions = self.init_decomposition(self.outer_polygon, bulk_reg, tol=self.mesh_step)
        self.group_positions[0] = np.mean(self.outer_polygon, axis=0)

        print("group positions ", self.group_positions)

        # extract fracture lines larger then the mesh step
        print("self.fr range ", self.fr_range)

        square_fr_range = [self.fr_range[0], np.min([self.fr_range[1], self.outer_polygon[1][0] - self.outer_polygon[0][0]])]
        print("square fr range ", square_fr_range)
        #self.fracture_lines = self.fractures.get_lines(self.fr_range)
        self.fracture_lines = self.fractures.get_lines(square_fr_range)
        print("len fracture lines ", len(self.fracture_lines))

        pd, fr_regions = self.add_fractures(pd, self.fracture_lines, eid=0)
        for reg in fr_regions:
            self.reg_to_group[reg.id] = 0
        self.decomp = pd

    def _make_mesh_window(self):
        print("make mesh window")
        import geometry_2d as geom

        self.skip_decomposition = False

        self.make_fracture_network()

        gmsh_executable = self.config_dict["gmsh_executable"]
        print("len self.regions ", len(self.regions))

        g2d = geom.Geometry2d("mesh_window_" + self.basename, self.regions)
        g2d.add_compoud(self.decomp)
        g2d.make_brep_geometry()
        print("self mesh step ", self.mesh_step)
        step_range = (self.mesh_step * 0.9, self.mesh_step * 1.1)
        g2d.call_gmsh(gmsh_executable, step_range)
        self.mesh = g2d.modify_mesh()

    def make_mesh(self):
        import geometry_2d as geom
        mesh_file = "mesh_{}.msh".format(self.basename)
        print("mesh file ", mesh_file)
        self.skip_decomposition = os.path.exists(mesh_file)

        print("self.skip decomposition ", self.skip_decomposition)

        #print("self regions ", self.regions)
        if not self.skip_decomposition:
            self.make_fracture_network()

            gmsh_executable = self.config_dict["gmsh_executable"]
            g2d = geom.Geometry2d("mesh_" + self.basename, self.regions)
            g2d.add_compoud(self.decomp)
            g2d.make_brep_geometry()
            step_range = (self.mesh_step * 0.9, self.mesh_step * 1.1)
            g2d.call_gmsh(gmsh_executable, step_range)
            self.mesh = g2d.modify_mesh()

        elif self.skip_decomposition and self.config_dict["mesh_window"]:
            self._make_mesh_window()
        else:
            self.make_fracture_network()
            self.mesh = gmsh_io.GmshIO()
            with open(mesh_file, "r") as f:
                self.mesh.read(f)

    def make_fields(self):
        """
        Calculate the conductivity and the cross-section fields, write into a GMSH file.

        :param cond_tensors_2d: Dictionary of the conductivities determined from a subscale
        calculation.
        :param cond_2d_samples: Array Nx2x2 of 2d tensor samples from own and other subsample problems.
        :return:
        """
        fracture_model = FractureModel(self.fractures, self.reg_to_fr, **self.config_dict['fracture_model'])
        elem_ids, cs_field, cond_tn_field = write_fields(self.mesh, self.basename, self.bulk_model, fracture_model)

        self._elem_ids = elem_ids
        self._cond_tn_field = cond_tn_field

    def bulk_field(self):
        assert self._elem_ids is not None

        points = []
        values00 = []
        values11 = []
        values01 = []
        for eid, tensor in zip(self._elem_ids, self._cond_tn_field):
            el_type, tags, node_ids = self.mesh.elements[eid]
            if len(node_ids) > 2:
                center = np.mean([np.array(self.mesh.nodes[nid]) for nid in node_ids], axis=0)
                points.append(center[0:2])
                # tensor is flatten 3x3
                values00.append(tensor[0])
                values01.append(tensor[1])
                values11.append(tensor[4])

        return np.array(points), np.array([values00, values01, values11])


    # def elementwise_mesh(self, coarse_mesh, mesh_step, bounding_polygon):
    #     import geometry_2d as geom
    #     from bgem.polygons.plot_polygons import plot_decomp_segments
    #
    #     mesh_file = "mesh_{}.msh".format(self.basename)
    #     if os.path.exists(mesh_file):
    #         # just initialize reg_to_group map
    #         self.skip_decomposition = True
    #
    #     self.mesh_step = mesh_step
    #     self.none_reg = self.add_region('none', dim=-1)
    #     # centers of macro elements
    #
    #     self.reg_to_group = {}  # bulk and fracture region id to coarse element id
    #     g2d = geom.Geometry2d("mesh_" + self.basename, self.regions, bounding_polygon)
    #     for eid, (tele, tags, nodes) in coarse_mesh.elements.items():
    #         # eid = 319
    #         # (tele, tags, nodes) = coarse_mesh.elements[eid]
    #         #print("Geometry for eid: ", eid)
    #         if tele != 2:
    #             continue
    #         prefix = "el_{:03d}_".format(eid)
    #         outer_polygon = np.array([coarse_mesh.nodes[nid][:2] for nid in nodes])
    #         # set mesh step to maximal height of the triangle
    #         area = np.linalg.norm(np.cross(outer_polygon[1] - outer_polygon[0], outer_polygon[2] - outer_polygon[0]))
    #
    #         self.mesh_step = min(self.mesh_step, area / np.linalg.norm(outer_polygon[1] - outer_polygon[0]))
    #         self.mesh_step = min(self.mesh_step, area / np.linalg.norm(outer_polygon[2] - outer_polygon[1]))
    #         self.mesh_step = min(self.mesh_step, area / np.linalg.norm(outer_polygon[0] - outer_polygon[2]))
    #
    #         self.group_positions[eid] = np.mean(outer_polygon, axis=0)
    #         #edge_sizes = np.linalg.norm(outer_polygon[:, :] - np.roll(outer_polygon, -1, axis=0), axis=1)
    #         #diam = np.max(edge_sizes)
    #
    #         bulk_reg = self.add_region(prefix + "bulk_2d_", dim=2, mesh_step=self.mesh_step)
    #         self.reg_to_group[bulk_reg.id] = eid
    #         # create regions
    #         # outer polygon
    #         normals = []
    #         shifts = []
    #         if eid==1353:
    #             print("break")
    #         pd, side_regions = self.init_decomposition(outer_polygon, bulk_reg, tol = self.mesh_step*0.8)
    #         for i_side, side_reg in enumerate(side_regions):
    #             side_reg.name = "." + prefix + side_reg.name[1:]
    #             side_reg.sub_reg.name = "." + prefix + side_reg.sub_reg.name[1:]
    #             self.reg_to_group[side_reg.id] = eid
    #             normals.append(side_reg.normal)
    #             shifts.append(side_reg.normal @ outer_polygon[i_side])
    #         self.side_regions.extend(side_regions)
    #
    #         # extract fracture lines larger then the mesh step
    #         fracture_lines = self.fractures.get_lines(self.fr_range)
    #         line_candidates = {}
    #         for i_fr, (p0, p1) in fracture_lines.items():
    #             for n, d in zip(normals, shifts):
    #                 sgn0 = n @ p0 - d > 0
    #                 sgn1 = n @ p1 - d > 0
    #                 #print(sgn0, sgn1)
    #                 if sgn0 and sgn1:
    #                     break
    #             else:
    #                 line_candidates[i_fr] = (p0, p1)
    #                 #print("Add ", i_fr)
    #             #plot_decomp_segments(pd, [p0, p1])
    #
    #
    #         pd, fr_regions = self.add_fractures(pd, line_candidates, eid)
    #         for reg in fr_regions:
    #             reg.name = prefix + reg.name
    #             self.reg_to_group[reg.id] = eid
    #         if not self.skip_decomposition:
    #             g2d.add_compoud(pd)
    #
    #     if self.skip_decomposition:
    #         self.mesh = gmsh_io.GmshIO()
    #         with open(mesh_file, "r") as f:
    #             self.mesh.read(f)
    #         return
    #
    #     g2d.make_brep_geometry()
    #     step_range = (self.mesh_step * 0.9, self.mesh_step *1.1)
    #     gmsh_executable = self.config_dict["gmsh_executable"]
    #     g2d.call_gmsh(gmsh_executable, step_range)
    #     self.mesh = g2d.modify_mesh()

    def run(self):
        outer_reg_names = []
        for reg in self.side_regions:
            outer_reg_names.append(reg.name)
            outer_reg_names.append(reg.sub_reg.name)
        self.thread = FlowThread(self.basename, outer_reg_names, self.config_dict)

        self.thread.start()
        return self.thread

    # def effective_tensor_from_balance(self, side_regions):
    #     """
    #     :param mesh: GmshIO mesh object.
    #     :param side_regions: List of side regions with the "normal" attribute.
    #     :return:
    #     """
    #     loads = self.pressure_loads
    #     with open(os.path.join(self.basename, "water_balance.yaml")) as f:
    #         balance = yaml.load(f, Loader=yaml.FullLoader)['data']
    #     flux_response = np.zeros_like(loads)
    #     reg_map = {}
    #     for reg in side_regions:
    #         reg_map[reg.name] = reg
    #         reg_map[reg.sub_reg.name] = reg
    #     for entry in balance:
    #         reg = reg_map.get(entry['region'], None)
    #         bc_influx = entry['data'][0]
    #         if reg is not None:
    #             flux_response[entry['time']] += reg.normal * bc_influx
    #     flux_response /= len(side_regions)
    #     #flux_response *= np.array([100, 1])[None, :]
    #
    #
    #     # least square fit for the symmetric conductivity tensor
    #     rhs = flux_response.flatten()
    #     # columns for the tensor values: C00, C01, C11
    #     pressure_matrix = np.zeros((len(rhs), 3))
    #     for i_load, (p0, p1) in enumerate(loads):
    #         i0 = 2 * i_load
    #         i1 = i0 + 1
    #         pressure_matrix[i0] = [p0, p1, 0]
    #         pressure_matrix[i1] = [0, p0, p1]
    #     C = np.linalg.lstsq(pressure_matrix, rhs, rcond=None)[0]
    #     cond_tn = np.array([[C[0], C[1]], [C[1], C[2]]])
    #     self.plot_effective_tensor(flux_response, cond_tn, label)
    #     #print(cond_tn)
    #     return cond_tn

    def element_volume(self, mesh, nodes):
        nodes = np.array([mesh.nodes[nid] for nid in  nodes])
        if len(nodes) == 1:
            return 0
        elif len(nodes) == 2:
            return np.linalg.norm(nodes[1] - nodes[0])
        elif len(nodes) == 3:
            return np.linalg.norm(np.cross(nodes[1] - nodes[0], nodes[2] - nodes[0]))
        else:
            assert False

    def effective_tensor_from_bulk(self):
        """
        :param bulk_regions: mapping reg_id -> tensor_group_id, groups of regions for which the tensor will be computed.
        :return: {group_id: conductivity_tensor} List of effective tensors.
        """
        bulk_regions = self.reg_to_group
        print("bulk regions ", bulk_regions)

        out_mesh = gmsh_io.GmshIO()
        with open(os.path.join(self.basename, "flow_fields.msh"), "r") as f:
            out_mesh.read(f)
        time_idx = 0
        time, field_cs = out_mesh.element_data['cross_section'][time_idx]
        ele_reg_vol = {eid: (tags[0] - 10000, self.element_volume(out_mesh, nodes))
                       for eid, (tele, tags, nodes) in out_mesh.elements.items()}

        assert len(field_cs) == len(ele_reg_vol)
        velocity_field = out_mesh.element_data['velocity_p0']

        print("velocity field ", velocity_field)

        loads = self.pressure_loads
        group_idx = {group_id: i_group for i_group, group_id in enumerate(set(bulk_regions.values()))}
        n_groups = len(group_idx)
        group_labels = n_groups * ['_']
        for reg_id, group_id in bulk_regions.items():
            i_group = group_idx[group_id]
            old_label = group_labels[i_group]
            new_label = self.regions[reg_id].name
            group_labels[i_group] = old_label if len(old_label) > len(new_label) else new_label

        n_directions = len(loads)
        flux_response = np.zeros((n_groups, n_directions, 2))
        area = np.zeros((n_groups, n_directions))
        print("Averaging velocities ...")
        for i_time, (time, velocity) in velocity_field.items():
            for eid, ele_vel in velocity.items():
                reg_id, vol = ele_reg_vol[eid]
                cs = field_cs[eid][0]
                volume = cs * vol
                i_group = group_idx[bulk_regions[reg_id]]
                flux_response[i_group, i_time, :] += -(volume * np.array(ele_vel[0:2]))
                area[i_group, i_time] += volume
        flux_response /= area[:, :, None]
        cond_tensors = {}
        print("Fitting tensors ...")
        for group_id, i_group in group_idx.items():
            flux = flux_response[i_group]
            # least square fit for the symmetric conductivity tensor
            rhs = flux.flatten()
            # columns for the tensor values: C00, C01, C11
            pressure_matrix = np.zeros((len(rhs), 3))
            for i_load, (p0, p1) in enumerate(loads):
                i0 = 2 * i_load
                i1 = i0 + 1
                pressure_matrix[i0] = [p0, p1, 0]
                pressure_matrix[i1] = [0, p0, p1]
            C = np.linalg.lstsq(pressure_matrix, rhs, rcond=None)[0]
            cond_tn = np.array([[C[0], C[1]], [C[1], C[2]]])
            if i_group < 10:
                 # if flux.shape[0] < 5:
                 print("Plot tensor for eid: ", group_id)
                 print("Fluxes: \n", flux)
                 print("pressures: \n", loads)
                 print("cond: \n", cond_tn)
                 self.plot_effective_tensor(flux, cond_tn, self.basename + "_" + group_labels[i_group])
                 #print(cond_tn)
            cond_tensors[group_id] = cond_tn
        self.cond_tensors = cond_tensors
        return cond_tensors

    def labeled_arrow(self, ax, start, end, label):
        """
        Labeled and properly scaled arrow.
        :param start: origin point, [x,y]
        :param end: tip point [x,y]
        :param label: string label, placed near the tip
        :return:
        """
        scale = np.linalg.norm(end - start)
        ax.arrow(*start, *end, width=0.003 * scale, head_length=0.1 * scale, head_width =0.05 * scale)
        if (end - start)[1] > 0:
            vert_align = 'bottom'
        else:
            vert_align = 'top'
        ax.annotate(label, end + 0.1*(end - start), va=vert_align)

    def plot_effective_tensor(self, fluxes, cond_tn, label):
        """
        Plot response fluxes for pressure gradients with angles [0, pi) measured from the X-azis counterclockwise.
        :param fluxes: np.array of 2d fluex vectors.
        :return:
        """

        import matplotlib.pyplot as plt

        e_val, e_vec = np.linalg.eigh(cond_tn)
        fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(10, 10))
        ax = axes
        # setting the axis limits in [left, bottom, width, height]
        #rect = [0.1, 0.1, 0.8, 0.8]
        ax.set_aspect('equal')
        #ax_polar = fig.add_axes(rect, polar=True, frameon=False)
        continuous_loads = np.array([(np.cos(t), np.sin(t)) for t in np.linspace(0, np.pi, 1000)])
        X, Y = cond_tn @ continuous_loads.T
        # print("Fluxes: ", fluxes)
        ax.scatter(X, Y, c='green', s=0.1)
        ax.scatter(-X, -Y, c='green', s=0.1)
        ax.scatter(fluxes[:, 0], fluxes[:, 1], c='red', s=30, marker='+')
        ax.scatter(-fluxes[:, 0], -fluxes[:, 1], c='red', s=30, marker='+')
        lim = max(max(X), max(fluxes[:, 0]), max(Y), max(fluxes[:, 1])) * 1.2
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        self.labeled_arrow(ax, [0,0], 0.9 * e_val[0] * e_vec[:, 0], "{:5.2g}".format(e_val[0])  )
        self.labeled_arrow(ax, [0,0], 0.9 * e_val[1] * e_vec[:, 1], "{:5.2g}".format(e_val[1])  )
        fig.suptitle("Conductivity tensor: {}".format(label))

        #ax_polar.grid(True)
        fig.savefig("cond_tn_{}.pdf".format(label))
        plt.close(fig)
        #plt.show()

    def summary(self):
        return dict(
            pos=[self.group_positions[eid].tolist() for eid in self.cond_tensors.keys()],
            cond_tn=[self.cond_tensors[eid].tolist() for eid in self.cond_tensors.keys()]
        )


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
        # sample_config attributes:
        # finer_level_path - Path to the file with microscale tensors from the previous level. Used for sampling conductivity.
        # config_path
        # do_coarse
        # h_coarse_step
        # h_fine_step
        # seedm
        # i_level
        print("sample config ", sample_config)

        self.__dict__.update(sample_config)

        print("self levels ", self.levels)

        self.i_level = 1
        self.finer_level_path = None  # It should be set if config param 'choose_from_finer_level' is True
        self.h_fine_step = self.levels[1]["step"]
        self.h_coarse_step = self.levels[0]["step"]

        self.do_coarse = False

        np.random.seed(self.seed)
        # with open(self.config_path, "r") as f:
        #     self.config_dict = yaml.load(f) # , Loader=yaml.FullLoader
        self.config_dict = sample_config

    def generate_fractures(self):
        geom = self.config_dict["geometry"]
        lx, ly = geom["fractures_box"]
        fr_size_range = geom["pow_law_size_range"]
        pow_law_exp_3d = geom["pow_law_size_exp"]
        pow_law_sample_range = geom["pow_law_sample_range"]
        n_frac_limit = geom["n_frac_limit"]
        print("n frac limit ", n_frac_limit)
        print("fr size range ", fr_size_range)
        p_32 = geom["p_32"]

        # generate fracture set
        fracture_box = [lx, ly, 0]
        area = lx * ly

        pop = fracture.Population(area, fracture.LineShape)
        pop.add_family("all",
                       fracture.FisherOrientation(0, 90, np.inf),
                       fracture.VonMisesOrientation(0, 0),
                       fracture.PowerLawSize.from_mean_area(pow_law_exp_3d - 1, fr_size_range, p_32, pow_law_exp_3d)
                       )

        if pow_law_sample_range:
            pop.set_sample_range(pow_law_sample_range)
        elif n_frac_limit:
            # pop.set_sample_range([None, np.min([lx, ly, self.config_dict["geometry"]["fr_max_size"]])],
            #                      sample_size=n_frac_limit)
            # pop.set_sample_range([None, max(lx, ly)],
            #                      sample_size=n_frac_limit)

            pop.set_sample_range(fr_size_range,
                                 sample_size=n_frac_limit)

        print("total mean size: ", pop.mean_size())
        print("size range:", pop.families[0].size.sample_range)

        print("pop families ", pop.families[0].size)


        pos_gen = fracture.UniformBoxPosition(fracture_box)
        fractures = pop.sample(pos_distr=pos_gen, keep_nonempty=True)

        print("fractures len ", len(fractures))

        fr_set = fracture.Fractures(fractures, fr_size_range[0] / 2)
        return fr_set

    def make_summary(self, done_list):
        results = {problem.basename: problem.summary() for problem in done_list}

        with open("summary.yaml", "w") as f:
            yaml.dump(results, f)

    def calculate(self):
        fractures = self.generate_fractures()
        print("fractures ", fractures)

        # fine problem
        fine_flow = FlowProblem.make_fine(self.i_level, (self.h_fine_step, self.config_dict["geometry"]["fr_max_size"]), fractures, self.finer_level_path,
                                          self.config_dict)
        fine_flow.make_mesh()
        fine_flow.make_fields()
        fine_flow.run()
        done = []
        # coarse problem
        # if self.do_coarse:
        #     coarse_ref = FlowProblem.make_microscale(self.i_level, (self.h_fine_step, self.h_coarse_step), fractures, fine_flow, self.config_dict)
        #     coarse_flow = FlowProblem.make_coarse(self.i_level, (self.h_coarse_step, np.inf), fractures, coarse_ref, self.config_dict)
        #
        #     # coarse mesh
        #     coarse_flow.make_fracture_network()
        #     coarse_flow.make_mesh()
        #
        #     # microscale mesh and run
        #     #coarse_ref.make_fracture_network()
        #     coarse_ref.elementwise_mesh(coarse_flow.mesh, self.h_fine_step,  coarse_flow.outer_polygon)
        #     coarse_ref.make_fields()
        #     coarse_ref.run().join()
        #     done.append(coarse_ref)
        #
        #     # coarse fields and run
        #     coarse_flow.make_fields()
        #     coarse_flow.run()
        #     coarse_flow.thread.join()
        #     done.append(coarse_flow)
        #     coarse_flow.effective_tensor_from_bulk()
        fine_flow.thread.join()
        done.append(fine_flow)
        fine_flow.effective_tensor_from_bulk()
        self.make_summary(done)

    # def mean_tensor(self, balance_dict, normals, regions):
    #     for n, reg in zip(normals, regions):
    #         balance_dict


def finished(start_time):
    sample_time = time.time() - start_time
    time.sleep(1)
    with open("FINISHED", "w") as f:
        f.write(f"done\n{sample_time}")


def compute_variogram(cond_field_xy, scalar_cond_log, dir):
    import matplotlib.pyplot as plt
    import itertools
    import skgstat as skg

    # Select pairs to sample various point distances
    #radius = 0.5 * np.linalg.norm(points.max_pt - points.min_pt, 2)
    n_samples = len(cond_field_xy)
    assert len(scalar_cond_log) == n_samples

    #all_combinations = list(itertools.combinations(cond_field_xy, 2))

    cond_field_xy_arr = np.squeeze(np.array(cond_field_xy))
    scalar_cond_log_arr = np.squeeze(np.array(scalar_cond_log))

    print(cond_field_xy_arr.shape)

    print(scalar_cond_log_arr.shape)

    #print("cond_field_values_arr[:,0] ", cond_field_values_arr[:, 0,0].shape)


    #V1 = skg.Variogram(cond_field_xy_arr, scalar_cond_log_arr, use_nugget = True, n_lags=100)
    #V1 = skg.Variogram(cond_field_xy_arr, scalar_cond_log_arr, n_lags=50, maxlag=0.8)
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    fig2, ax2 = plt.subplots(1, 1, figsize=(8, 6))
    V1 = skg.Variogram(cond_field_xy_arr, scalar_cond_log_arr)
    v1_figure = V1.plot(axes=ax, hist=True)
    print("v1 figure ", v1.figure)
    fig.savefig(os.path.join(dir, "variogram.pdf"))
    #V1.scattergram()
    V1.distance_difference_plot(ax=ax2)
    fig2.savefig(os.path.join(dir, "dist_differences.pdf"))




    # print("V1 coordinates ", V1.coordinates)
    #
    # print("len(lag groups)", len(V1.lag_groups()))
    # print("lag groups", V1.lag_groups())
    # print("V1.data()", V1.data())
    # print("V1.bins", V1.bins)
    # print("V1.diff() ", V1._diff)
    # print("V.metric_space ", V1.metric_space.dist_metric)
    #
    #
    # print("mode V1 distance matrix ", stats.mode(V1.distance_matrix, axis=None))
    #
    # distances = distance.cdist(cond_field_xy_arr, cond_field_xy_arr, 'euclidean')
    # print("distances ", distances)
    # print("len(distances.flatten) ", len(distances.flatten()))
    # print("len distances ", len(distances))
    #
    # print("len V1.data ", len(V1.data()))
    #
    # for lag_class in V1.lag_classes():
    #     print("lag class ", lag_class)
    # exit()

    ## Directional variogram
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    Vnorth = skg.DirectionalVariogram(cond_field_xy_arr, scalar_cond_log_arr, azimuth=90, tolerance=90)
    Vnorth.plot(axes=ax,show=True)
    Vnorth.distance_difference_plot()
    fig.savefig(os.path.join(dir,"variogram_vnorth.pdf"))

    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    Veast = skg.DirectionalVariogram(cond_field_xy_arr, scalar_cond_log_arr, azimuth=0, tolerance=90)
    Veast.plot(axes=ax, show=True)
    Veast.distance_difference_plot()
    fig.savefig(os.path.join(dir,"variogram_veast.pdf"))

    # V2 = skg.DirectionalVariogram(cond_field_xy_arr, scalar_cond_log_arr, azimuth=45, tolerance=30)
    # V2.plot(show=True)

    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    ax.plot(Vnorth.bins, Vnorth.experimental, '.--r', label='North-South')
    ax.plot(Veast.bins, Veast.experimental, '.--b', label='East-West')
    #ax.plot(V2.bins, V2.experimental, '.--g')
    ax.set_xlabel('lag [m]')
    ax.set_ylabel('semi-variance (matheron)')
    plt.legend(loc='upper left')
    fig.show()
    fig.savefig(os.path.join(dir, "directional_semi_variance.pdf"))

    exit()


def process_results(sample_dict):
    seeds = sample_dict.get("seeds", [])
    if len(seeds) == 0:
        seeds = [12345]

    seeds_cond_field_xy = []
    seeds_scalar_cond_log = []
    for seed in seeds:
        i = 0
        cond_field_xy = []
        scalar_cond_log = []

        subdir = "n_10_s_300_300"

        dir = "seed_{}/{}".format(seed, subdir)

        if not os.path.exists(dir):
            print("dir not exists: {}".format(dir))
            continue

        while True:
            dir_name = os.path.join(dir, "fine_{}".format(i))
            print("cwd ", os.getcwd())
            print("dir name ", dir_name)

            if os.path.exists(dir_name):
                print("dir exists ")
                with open(os.path.join(dir_name, "summary.yaml"), "r") as f:
                    summary_dict = yaml.load(f, Loader=yaml.FullLoader)
                print("   ...store")
                fine_cond_tn = np.array(summary_dict['fine']['cond_tn'])
                e_val, e_vec = np.linalg.eigh(fine_cond_tn)
                cond_field_xy.append(np.array(summary_dict['fine']['pos']))
                scalar_cond_log.append(np.log(np.mean(e_val)))
            else:
                break
            i+=1

        seeds_cond_field_xy.append(cond_field_xy)
        seeds_scalar_cond_log.append(scalar_cond_log)

    print("len(seeds_cond_field_xy) ", len(seeds_cond_field_xy))

    if len(seeds_cond_field_xy) > 0:  # List of samples, every sample have conductivity tensor for every coarse mesh element.
        for idx, (s_cond_field_xy, s_scalar_cond_log) in enumerate(zip(seeds_cond_field_xy, seeds_scalar_cond_log)):
            # self.precompute()
            # self.test_homogenity()
            # self.test_homogenity_xy()
            # self.test_isotropy()
            # self.test_isotropy_alt()
            print("s cond field xy ", s_cond_field_xy)
            print("s scalar cond log ", s_scalar_cond_log)

            compute_variogram(s_cond_field_xy, s_scalar_cond_log, dir="seed_{}/{}".format(seeds[idx], subdir))

            # self.eigenvals_correlation()
            # self.calculate_field_parameters()

    compute_anova(seeds_scalar_cond_log)

    print("seeds cond field xy ", seeds_cond_field_xy)
    print("seeds scalar cond log ", seeds_scalar_cond_log)


def compute_anova(seeds_scalar_cond_log):
    import scipy.stats as stats
    print("compute anova")
    # stats f_oneway functions takes the groups as input and returns ANOVA F and p value
    print("seeds_scalar_cond_log shape ", np.array(seeds_scalar_cond_log).shape)

    seeds_scalar_cond_log_arr = np.array(seeds_scalar_cond_log)
    print("seeds_scalar_cond_log_arr", seeds_scalar_cond_log_arr.shape)
    groups = []
    for i in range(seeds_scalar_cond_log_arr.shape[1]):
        groups.append(seeds_scalar_cond_log_arr[:, i])

    print("groups shape ", np.array(groups).shape)
    print(groups)

    fvalue, pvalue = stats.f_oneway(*groups) # (group 1 samples, group 2 sampels, ...)
    print(fvalue, pvalue)

    significant_level = 0.05
    if pvalue > significant_level:
        print("Cannot reject the null hypothesis that the population means are all equal")
    else:
        print("The null hypothesis is rejected - not all of population means are equal.")

    exit()


def subdomains(work_dir, sample_dict, command="run"):

    if command == "process":
        process_results(work_dir, sample_dict)
    else:
        run_samples(work_dir, sample_dict)


def run_samples(work_dir, sample_dict):
    domain_box = sample_dict["geometry"]["domain_box"]
    subdomain_box = sample_dict["geometry"]["subdomain_box"]
    subdomain_overlap = np.array([0, 0])  # np.array([50, 50])

    print("domain box ", domain_box)
    print("subdomain box ", subdomain_box)

    lx, ly = domain_box

    bottom_left_corner = [-lx / 2, -ly / 2]
    bottom_right_corner = [+lx / 2, -ly / 2]
    top_right_corner = [+lx / 2, +ly / 2]
    top_left_corner = [-lx / 2, +ly / 2]
    #               bottom-left corner  bottom-right corner  top-right corner    top-left corner
    complete_polygon = [bottom_left_corner, bottom_right_corner, top_right_corner, top_left_corner]

    plt.scatter(*zip(*complete_polygon))

    n_subdomains = sample_dict["geometry"].get("n_subdomains", 4)

    domain_box = sample_dict["geometry"]["domain_box"]
    subdomain_box = sample_dict["geometry"]["subdomain_box"]
    lx, ly = domain_box

    k = 0
    for i in range(n_subdomains):
        center_x = subdomain_box[0] / 2 + (lx - subdomain_box[0]) / (n_subdomains - 1) * i - lx / 2
        print("center_x ", center_x)
        for j in range(n_subdomains):
            k += 1
            print("i: {}, j:{}, k:{}".format(i, j, k))
            # if k != 24:
            #     continue
            if k != 6: #@TODO: zkusit jeste zmensit velikost elementu
                continue
            center_y = subdomain_box[1] / 2 + (lx - subdomain_box[1]) / (n_subdomains - 1) * j - lx / 2

            bl_corner = [center_x - subdomain_box[0] / 2, center_y - subdomain_box[1] / 2]
            br_corner = [center_x + subdomain_box[0] / 2, center_y - subdomain_box[1] / 2]
            tl_corner = [center_x - subdomain_box[0] / 2, center_y + subdomain_box[1] / 2]
            tr_corner = [center_x + subdomain_box[0] / 2, center_y + subdomain_box[1] / 2]

            outer_polygon = [copy.deepcopy(bl_corner), copy.deepcopy(br_corner), copy.deepcopy(tr_corner),
                             copy.deepcopy(tl_corner)]
            print("outer polygon ", outer_polygon)


            plt.scatter(*zip(*outer_polygon))

            sample_dict["geometry"]["outer_polygon"] = outer_polygon

            sample_dict["work_dir"] = work_dir

            bs = BothSample(sample_dict)
            bs.calculate()

            dir_name = "fine_{}".format(k)
            sample_dict["dir_name"] = dir_name
            try:
                import shutil
                shutil.move("fine", dir_name)
                shutil.move("summary.yaml", dir_name)
                shutil.move("mesh_fine.msh", dir_name)
                shutil.move("mesh_fine.brep", dir_name)
                shutil.move("mesh_fine.tmp.geo", dir_name)
                shutil.move("mesh_fine.tmp.msh", dir_name)
                shutil.rmtree("fine")
            except:
                pass

    plt.show()

    # bottom_left_corner = [-lx / 2, -ly / 2]
    # bottom_right_corner = [+lx / 2, -ly / 2]
    # top_right_corner = [+lx / 2, +ly / 2]
    # top_left_corner = [-lx / 2, +ly / 2]
    # #               bottom-left corner  bottom-right corner  top-right corner    top-left corner
    # complete_polygon = [bottom_left_corner, bottom_right_corner, top_right_corner, top_left_corner]
    # # bs = BothSample(sample_dict)
    # # bs.calculate()
    #
    # plt.scatter(*zip(*complete_polygon))
    #
    # tl_corner = top_left_corner
    #
    # i = 0
    # while True:
    #     bl_corner = tl_corner - np.array([0, subdomain_box[1]])
    #     bl_corner = np.maximum(bl_corner,
    #                            [bl_corner[0], bottom_right_corner[1]])  # make sure it lies in the complete polygon
    #
    #     br_corner = bl_corner + np.array([subdomain_box[0], 0])
    #     br_corner = np.minimum(br_corner, [bottom_right_corner[0], br_corner[1]])
    #
    #     tr_corner = tl_corner + np.array([subdomain_box[0], 0])
    #     tr_corner = np.minimum(tr_corner, [bottom_right_corner[0], tr_corner[1]])
    #
    #     outer_polygon = [copy.deepcopy(bl_corner), copy.deepcopy(br_corner), copy.deepcopy(tr_corner),
    #                      copy.deepcopy(tl_corner)]
    #     print("outer polygon ", outer_polygon)
    #
    #     sample_dict["geometry"]["outer_polygon"] = outer_polygon
    #
    #     # plt.scatter(*zip(*outer_polygon))
    #
    #     # Set new top-left corner including overlap on x axes
    #     # Do not move corner if at the end of x axis
    #     if tr_corner[0] - tl_corner[0] == subdomain_box[0]:
    #         tl_corner = tr_corner - np.array([subdomain_overlap[0], 0])
    #     else:
    #         tl_corner = tr_corner
    #
    #     # top-left corner is on the left-hand edge
    #     if tl_corner[0] >= bottom_right_corner[0]:
    #         tl_corner[0] = bottom_left_corner[0]
    #         tl_corner[1] = tl_corner[1] - subdomain_box[1]
    #         # y-axis overlap
    #         tl_corner += np.array([0, subdomain_overlap[1] if tl_corner[1] < top_left_corner[1] else 0])
    #
    #     # bs = BothSample(sample_dict)
    #     # bs.calculate()
    #
    #     dir_name = "fine_{}".format(i)
    #     sample_dict["dir_name"] = dir_name
    #     try:
    #         os.remove("mesh_fine.msh")
    #         os.remove("mesh_fine.brep")
    #         os.remove("mesh_fine.tmp.geo")
    #         os.remove("mesh_fine.tmp.msh")
    #         import shutil
    #         shutil.move("fine", dir_name)
    #         shutil.copy("summary.yaml", dir_name)
    #         shutil.rmtree("fine")
    #     except:
    #         pass
    #
    #
    #     # top-left corner is on the bottom edge
    #     if tl_corner[1] <= bottom_left_corner[1]:
    #         break
    #
    #     i += 1
    #
    # plt.show()



if __name__ == "__main__":
    start_time = time.time()
    atexit.register(finished, start_time)
    work_dir = sys.argv[1]
    sample_config = sys.argv[2]
    #try:
    with open(sample_config, "r") as f:
        sample_dict = yaml.load(f, Loader=yaml.FullLoader)
    #except Exception as e:
    #    print("cwd: ", os.getcwd(), "sample config: ", sample_config)

    subdomains(work_dir, sample_dict, command="run")
    #subdomains(sample_dict, command="process")

