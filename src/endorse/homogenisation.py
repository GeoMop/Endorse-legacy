
"""
Implement static homogenisation.
Setup an object collection subproblems, covering some domain.
Do subdomains calculations for various boundary conditions. These sampling could be done iteratively.
For given set of evaluation points, evaluate the properties in these points:
1. find appropriate subdomain with max. intersection with the negbourhood
   This step provides topological ingormation that is independent of the evaluatied quantities.
2. Compute ekvivalent properties in the points:
   scalar - just volume avarage with kernel weight
   vector - just volume avarage with kernel weight
   tensor - implicitely by 2 vector quantities
   tensor 4 - implicitely by 2 tensor quantities (could be outer product of vectors)

   symmetric tensor from outer product: a2+b1 , a3+c1, b3+c2
"""
import logging
from typing import *
import os
import numpy as np
from bgem.gmsh.gmsh_io import GmshIO
import bih
import attrs
import copy

from .mesh_class import Mesh, Element
from . import common
from .common import dotdict, memoize, File, report


@attrs.define
class Subdomain:
    mesh: Mesh
    el_indices: List[int]

    @staticmethod
    def for_element(micro_mesh: Mesh, macro_el: Element):
        """
        Select elements from the micro mesh interacting with a sphere
        approximating the macro element `id_el`.
        """
        center = macro_el.barycenter()
        distances = np.linalg.norm(macro_el.vertices() - center[None,:], axis=1)
        r = np.mean(distances)
        #logging.info(f"{center}, {distances}, {r}")
        return Subdomain.for_sphere(micro_mesh, center, r)

    @staticmethod
    def interact_sphere(center, r, nodes:np.array):
        """
        If element has any node in the subdomain.
        """
        nc = nodes[:, :] - center[None, :]
        indicate = np.sum(nc * nc, axis=1) < r*r    # shape n - nodes
        return np.any(indicate)


    @staticmethod
    def for_sphere(mesh, center, r) -> 'Subdomain':
        center = np.atleast_1d(center)
        aabb = bih.AABB([center - r, center + r])
        candidates = mesh.candidate_indices(aabb)
        assert candidates, f"Sphere AABB: {aabb} out of mesh AABB: {mesh.bih.aabb()}"
        subdomain_indices = [ie for ie in candidates
                         if Subdomain.interact_sphere(center, r, mesh.elements[ie].vertices())]
        #logging.info(f"Subdomain candidates: {len(candidates)}, elements: {len(subdomain_indices)}")
        assert subdomain_indices
        return Subdomain(mesh, subdomain_indices)

    def average(self, element_vec_data: np.array):
        """
        :param element_vec_data: Values on all self.mesh elements in the same ordering (checked in get_statice
        :return:
        """
        if len(element_vec_data.shape) == 1:
            element_vec_data = element_vec_data[:, None]
        assert len(self.mesh.elements) == element_vec_data.shape[0]
        volumes = self.mesh.el_volumes[self.el_indices]
        sub_domain_vector = element_vec_data[self.el_indices, :]
        avg = np.sum(sub_domain_vector[:, :] * volumes[:,None], axis=0) / np.sum(volumes)
        if avg.shape[0] == 1:
            return avg[0]
        else:
            return avg

# def micro_response(subdomains):
#     mesh = GmshIO("output/flow_fields.msh")
#     tree, el_ids = mesh_build_bih(mesh)
#     for x,y,z,r in subdomains:
#         center = np.array([x,y,z])
#         box = bih.AABB([center - r, center + r])
#         candidates = tree.find_box(box)
#         subdomain_els = [subdomain_interract for iel in candidates:


def make_subdomains_old(cfg, subdomains_params):
    fractures = []
    i_pos = 0
    mesh_file = "transport_micro.msh"
    fine_micro_mesh(cfg.geometry, fractures, i_pos, mesh_file)
    fine_mesh = Mesh(GmshIO(mesh_file))
    subdomains = [Subdomain.for_sphere(fine_mesh, center, r)
                  for center, r in subdomains_params]
    return subdomains

@memoize
def subdomains_mesh(subdomains: List[Subdomain], output_name="fine_with_subdomains.msh"):
    """
    Fine mesh with duplicated elements for the subdomains,
    markeg by regions "sub_<isubdomain>"
    """
    mesh = subdomains[0].mesh
    output_mesh = copy.copy(mesh.gmsh_io)
    el_id = max(mesh.el_ids)
    for isub, sub in enumerate(subdomains):
        reg_name = f"sub_{isub}"
        tag = 1000 + isub
        output_mesh.physical[reg_name] = (tag, 3)
        els = {}
        for ie in sub.el_indices:
            type, tags, inodes = output_mesh.elements[mesh.el_ids[ie]]
            physical_tag, entity_tag = tags
            el_id = el_id + 1
            els[el_id] = (type, (tag, tag), inodes)
        output_mesh.elements.update(els)
    output_mesh.write(output_name)



class Homogenisation:
    """
    Support for (implicit) homogenisation. General algorithm:
    1. create a micro problem mesh
    2. for this micro mesh create subdomains (subsets of micro elements)
       could possibly be weighted by suitable kernel see Subdomain class
    3. create Homogenisation for the micro mesh and defined subdomains.
    4. Run micro problem for some set of boundary or initial conditions,
       use Homogenisation methods to compute Subdomain averages of given scalar or vector quantities.
       Not stored but returned as numpy arrays.
    5. Use Homogenisation method to calculate ekvivalent tensor field for pair of vector fields.
    6. scalar quantities, meight be homogenised without running the micro problem, but we support them anyway.

    TODO: improve design to naturaly satisfy usage constrains:
      - common mesh for subdomains
      - gmsh_io for compatible mesh as the input !!
      - compatible response and load samples (comming from same solution)
      - rather functional design, homogenisation just as the subdomain averaging tool,
        for every sample apply custom averaging to get pairs fo load and response averages (like current impl)
        Finally apply Homogenisation method to derive fled of ekvivalent tensors.
        ... Homogenisation does not store the data so it does not modify.
    """
    def __init__(self, subdomains: List[Subdomain]):
        self.subdomains = subdomains
        self.micro_mesh = self.subdomains[0].mesh
        # We assume common mesh for all subdomains, TODO: better design to have this automaticaly satisfied.
        assert all([sub.mesh is self.micro_mesh for sub in self.subdomains])

    def average(self, element_wise_field):
        return np.array([sub.average(element_wise_field) for sub in self.subdomains])


    @staticmethod
    def equivalent_tensor_3d(loads, responses):
        # tensor pos. def.  <=> load @ response > 0
        # ... we possibly modify responses to satisfy
        unit_loads = loads / np.linalg.norm(loads, axis=1)[:, None]
        load_components = np.sum(responses * unit_loads, axis=1)
        responses_fixed = responses + (np.maximum(0, load_components) - load_components)[:, None] * unit_loads
        # from LS problem for 6 unknowns in Voigt notation: X, YY, ZZ, YZ, XZ, XY
        # the matrix has three blocks for Vx, Vy, Vz component of the responses
        # each block has different sparsity pattern
        n_loads = loads.shape[0]
        zeros = np.zeros(n_loads)
        ls_mat_vx = np.stack([loads[:, 0], zeros, zeros, zeros, loads[:, 2], loads[:, 1]], axis=1)
        rhs_vx = responses_fixed[:, 0]
        ls_mat_vy = np.stack([zeros, loads[:, 1], zeros, loads[:, 2], zeros, loads[:, 0]], axis=1)
        rhs_vy = responses_fixed[:, 1]
        ls_mat_vz = np.stack([zeros, zeros, loads[:, 2], loads[:, 1], loads[:, 0], zeros], axis=1)
        rhs_vz = responses_fixed[:, 2]
        ls_mat = np.concatenate([ls_mat_vx, ls_mat_vy, ls_mat_vz], axis=0)
        rhs = np.concatenate([rhs_vx, rhs_vy, rhs_vz], axis=0)
        assert ls_mat.shape == (3 * n_loads, 6)
        assert rhs.shape == (3 * n_loads,)
        result = np.linalg.lstsq(ls_mat, rhs, rcond=None)
        cond_tn_voigt, residuals, rank, singulars = result
        condition_number = singulars[0] / singulars[-1]
        if condition_number > 1e3:
            logging.warning(f"Badly conditioned inversion. Residual: {residuals}, max/min sing. : {condition_number}")
        return cond_tn_voigt

    #@report
    def equivalent_tensor_field(self, load_field,  response_field):
        load_field = np.array(load_field)
        response_field = np.array(response_field)
        assert load_field.shape == response_field.shape
        assert load_field.shape[1] == len(self.subdomains)
        assert load_field.shape[2] == 3
        tensors = np.empty((len(self.subdomains), 6))
        for isub, sub in enumerate(self.subdomains):
            loads = load_field[:, isub, :]
            responses = response_field[:, isub, :]
            tensors[isub, :] = self.equivalent_tensor_3d(loads, responses)
        return tensors








#
# def eval_conductivity_field(cfg_micro_conductivity, eval_points):
#     fine_mesh, max_radius, output_file):
#     c_min, c_max = cfg_micro_conductivity.range
#     x = bary_coords[:, 0]
#     y = bary_coords[:, 1]
#     z = bary_coords[:, 2]
#
#
#     horizontal_dist = y*y + (z*z*9)
#     vertical_dist = y*y*9 + z*z
#     r_sq = max_radius * max_radius
#     cond = c_min + c_max * np.exp(-horizontal_dist/r_sq) + c_max * np.exp(-vertical_dist/r_sq)
#
# def fine_conductivity_field(cfg, fine_mesh, output_file):
#     bary_coords = np.array([el.barycenter() for el in fine_mesh.elements])
#     conductivity = eval_conductivity_field((cfg.fine_conductivity, bary_coords)
#     return fine_mesh.write_fields(output_file, dict(conductivity=conductivity))

# TODO:
# - review main structure from macro transport
# - introduce new main test
# - new subtests
#