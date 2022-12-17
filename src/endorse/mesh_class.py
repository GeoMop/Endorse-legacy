from typing import Dict, Tuple, List

import attrs
import bih
import numpy as np
from numba import njit

from bgem.gmsh.gmsh_io import GmshIO

from endorse.common import File, memoize, report


#@njit
def element_vertices(all_nodes: np.array, node_indices: np.array):
    return all_nodes[node_indices[:], :]


#@njit
def element_loc_mat(all_nodes: np.array, node_indices: List[int]):
    n = element_vertices(all_nodes, node_indices)
    return (n[1:, :] - n[0]).T


#@njit
def element_compute_volume(all_nodes: np.array, node_indices: List[int]):
    return np.linalg.det(element_loc_mat(all_nodes, node_indices)) / 6


@attrs.define
class Element:
    mesh: 'Mesh'
    type: int
    tags: Tuple[int, int]
    node_indices: List[int]

    def vertices(self):
        return element_vertices(self.mesh.nodes, np.array(self.node_indices, dtype=int))

    def loc_mat(self):
        return element_loc_mat(self.mesh.nodes, self.node_indices)

    def volume(self):
        return element_compute_volume(self.mesh.nodes, self.node_indices)

    def barycenter(self):
        return np.mean(self.vertices(), axis=0)

    def gmsh_tuple(self, node_map):
        node_ids = [node_map[inode] for inode in self.node_indices]
        return (self.type, self.tags, node_ids)



@memoize
def _load_mesh(mesh_file: File):
    # !! can not memoize static and class methods (have no name)
    return Mesh(GmshIO(mesh_file.path), file = mesh_file)


@report
#@njit
def mesh_compute_el_volumes(nodes:np.array, node_indices :np.array) -> np.array:
    return np.array([element_compute_volume(nodes, ni) for ni in node_indices])


class Mesh:

    @staticmethod
    def load_mesh(mesh_file: File) -> 'Mesh':
        return _load_mesh(mesh_file)

    @staticmethod
    def empty(mesh_path) -> 'Mesh':
        return Mesh(GmshIO(), mesh_path)

    def __init__(self, gmsh_io: GmshIO, file):

        self.gmsh_io : GmshIO = gmsh_io
        # TODO: remove relation to file
        # rather use a sort of generic wrapper around loadable objects
        # in order to relay on the underlaing files for the caching
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
        list_box = box.tolist()
        return self.bih.find_box(bih.AABB(list_box))

    # def el_volume(self, id):
    #     return self.elements[self.el_indices[id]].volume()

    @property
    @report
    def el_volumes(self):
        if self._el_volumes is None:
            node_indices = np.array([e.node_indices for e in self.elements], dtype=int)
            print(f"Compute el volumes: {self.nodes.shape}, {node_indices.shape}")
            self._el_volumes = mesh_compute_el_volumes(self.nodes, node_indices)
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

    def submesh(self, elements, file_path):
        gmesh = GmshIO()
        active_nodes = np.full( (len(self.nodes),), False)
        for iel in elements:
            el = self.elements[iel]
            active_nodes[el.node_indices] = True
        sub_nodes = self.nodes[active_nodes]
        new_for_old_nodes = np.zeros((len(self.nodes),), dtype=int)
        new_for_old_nodes[active_nodes] = np.arange(1,len(sub_nodes)+1, dtype=int)
        gmesh.nodes = {(nidx+1):node for nidx, node in enumerate(sub_nodes)}
        gmesh.elements = {(eidx+100): self.elements[iel].gmsh_tuple(node_map=new_for_old_nodes) for eidx, iel in enumerate(elements)}
        #print(gmesh.elements)
        gmesh.physical = self.gmsh_io.physical
        #gmesh.write(file_path)
        gmesh.normalize()
        return Mesh(gmesh, "")

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


    def write_fields(self, file_name:str, fields: Dict[str, np.array]=None) -> File:
        self.gmsh_io.write(file_name, format="msh2")
        if fields is not None:
            self.gmsh_io.write_fields(file_name, self.el_ids, fields)
        return File(file_name)

