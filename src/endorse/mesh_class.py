from typing import Dict, Tuple, List

import attrs
import bih
import numpy as np
from bgem.gmsh.gmsh_io import GmshIO

from endorse.common import File, memoize


@memoize
def _load_mesh(mesh_file: File):
    # !! can not memoize static and class methods (have no name)
    return Mesh(GmshIO(mesh_file.path), file = mesh_file)


class Mesh:

    @staticmethod
    def load_mesh(mesh_file: File) -> 'Mesh':
        return _load_mesh(mesh_file)

    def __init__(self, gmsh_io: GmshIO, file):

        self.gmsh_io : GmshIO = gmsh_io
        self.file : File = file
        self.reinit()


    def reinit(self):
        # bounding interval hierarchy for the mesh elements
        # numbers elements from 0 as they are added
        self.node_ids = []
        self.node_indices = {}
        self.nodes = np.empty((len(self.gmsh_io.nodes), 3))
        for i, (nid, node) in enumerate(self.gmsh_io.nodes.items()):
            self.node_indices[nid] = i
            self.node_ids.append(nid)
            self.nodes[i, :] = node

        self.el_ids = []
        self.el_indices = {}
        self.elements = []
        for i, (eid, el) in enumerate(self.gmsh_io.elements.items()):
            type, tags, node_ids = el
            element = Element(self, type, tags, [self.node_indices[nid] for nid in node_ids])
            self.el_indices[eid] = i
            self.el_ids.append(eid)
            self.elements.append(element)

        # _boxes: List[bih.AABB]
        self.bih: bih.BIH = self._build_bih()

        self._el_volumes:np.array = None
        self._el_barycenters:np.array =  None

    def __getstate__(self):
        return (self.gmsh_io, self.file)

    def __setstate__(self, args):
        self.gmsh_io, self.file = args
        self.reinit()

    def _build_bih(self):
        el_boxes = []
        for el in self.elements:
            node_coords = el.vertices()
            box = bih.AABB(node_coords)
            el_boxes.append(box)
        _bih = bih.BIH()
        _bih.add_boxes(el_boxes)
        _bih.construct()
        return _bih



    def candidate_indices(self, box):
        return self.bih.find_box(box)

    # def el_volume(self, id):
    #     return self.elements[self.el_indices[id]].volume()

    @property
    def el_volumes(self):
        if self._el_volumes is None:
            self._el_volumes = np.array([e.volume() for e in self.elements])
        return self._el_volumes

    def el_barycenters(self):
        if self._el_barycenters is None:
            self._el_barycenters = np.array([e.barycenter() for e in self.elements])
        return self._el_barycenters.T

    # def el_loc_mat(self, id):
    #     return self.elements[self.el_indices[id]].loc_mat()

    # def el_barycenter(self, id):
    #     return self.elements[self.el_indices[id]].barycenter()

    # def el_nodes(self, id):
    #     return self.elements[self.el_indices[id]].vertices()


    def get_static_p0_values(self, field_name:str):
        field_dict = self.gmsh_io.element_data[field_name]
        assert len(field_dict) == 1
        values = field_dict[0].values
        value_ids = field_dict[0].tags
        value_to_el_idx = [self.el_indices[iv] for iv in value_ids]
        values_mesh = np.empty_like(values)
        values_mesh[value_to_el_idx[:]] = values
        return values_mesh

    def get_static_p1_values(self, field_name:str):
        field_dict = self.gmsh_io.node_data[field_name]
        assert len(field_dict) == 1
        values = field_dict[0].values
        value_ids = field_dict[0].tags
        value_to_node_idx = [self.node_indices[iv] for iv in value_ids]
        values_mesh = np.empty_like(values)
        values_mesh[value_to_node_idx[:]] = values
        return values_mesh

    def write_fields(self, file_name:str, fields: Dict[str, np.array]) -> None:
        self.gmsh_io.write(file_name, format="msh2")
        self.gmsh_io.write_fields(file_name, self.el_ids, fields)
        return File(file_name)


@attrs.define
class Element:
    mesh: 'Mesh'
    type: int
    tags: Tuple[int, int]
    node_indices: List[int]

    def vertices(self):
        return self.mesh.nodes[self.node_indices, :]

    def loc_mat(self):
        n = self.vertices()
        return (n[1:, :] - n[0]).T

    def volume(self):
        return np.linalg.det(self.loc_mat()) / 6

    def barycenter(self):
        return np.mean(self.vertices(), axis=0)