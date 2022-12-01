import logging
import os
from typing import *

import numpy as np

from . import common
from .common import dotdict, memoize, File, call_flow, workdir, report
from .mesh import repository_mesh as repo_mesh, mesh_tools
from .homogenisation import  subdomains_mesh, Homogenisation, Subdomain, MacroSphere, make_subproblems, Subproblems
from .mesh.repository_mesh import one_borehole
from .mesh_class import Mesh
from . import apply_fields
from . import plots
from bgem.stochastic.fracture import Fracture

def input_files(cfg_tr_full):
    return [
        cfg_tr_full.piezo_head_input_file,
        cfg_tr_full.conc_flux_file
    ]

def fullscale_transport(cfg, source_params, seed):
    """
    1. apply conouctivity to given mesh:
       - on borehole neighbourhood, select elements
       - calculate barycenters
       - apply conductivity
       - write the field
    2. substitute source term space distribution
    3. return necessary files
    """

    cfg_fine = cfg.transport_fullscale
    full_mesh_file, fractures = fullscale_transport_mesh(cfg, cfg_fine.mesh, seed)

    full_mesh = Mesh.load_mesh(full_mesh_file, heal_tol = 1e-4)
    el_to_fr = fracture_map(full_mesh, fractures)
    #mesh_modified_file = full_mesh.write_fields("mesh_modified.msh2")
    #mesh_modified = Mesh.load_mesh(mesh_modified_file)

    input_fields_file = compute_fields(cfg, full_mesh, el_to_fr)
    large_model = File(os.path.basename(cfg_fine.piezo_head_input_file))
    conc_flux = File(os.path.basename(cfg_fine.conc_flux_file))
    params = cfg_fine.copy()
    new_params = dict(
        mesh_file=input_fields_file,
        piezo_head_input_file=large_model,
        conc_flux_file=conc_flux,
        input_fields_file = input_fields_file
    )
    params.update(new_params)
    params.update(set_source_limits(cfg))
    template = os.path.join(common.flow123d_inputs_path, cfg_fine.input_template)
    common.call_flow(cfg.flow_env, template, params)

def fracture_map(mesh, fractures) -> Dict[int, Fracture]:
    """
    - join all fracture regions into single "fractures" region
    - return dictionary mapping element idx to fracture
    :param mesh:
    :param fractures:
    :return:
    """
    own_name_to_id = {fr.region.name: fr.region.id for fr in fractures}
    own_to_gmsh_id = {own_name_to_id[name]: gmsh_id for name,(gmsh_id, dim) in mesh.gmsh_io.physical.items() if name in own_name_to_id}

    new_reg_id = max( [gmsh_id for gmsh_id, dim in mesh.gmsh_io.physical.values()] ) + 1
    new_reg_map = {(own_to_gmsh_id[fr.region.id], 2): (new_reg_id, 2, "fractures")  for fr in fractures}

    # if do_heal:
    #     hm.heal_mesh(gamma_tol=0.01)
    #     hm.move_all(geom_dict["shift_vec"])
    #     elm_to_orig_reg = hm.map_regions(new_reg_map)
    #     hm.stats_to_yaml(mesh_name + "_heal_stats.yaml")
    #     assert hm.healed_mesh_name == mesh_healed
    #     hm.write()
    # else:
    iel_to_orig_reg = mesh.map_regions(new_reg_map)

    reg_to_fr = {own_to_gmsh_id[fr.region.id]: fr for fr in fractures}
    elm_to_fr = {el_idx: reg_to_fr[reg_id] for el_idx, (reg_id, dim) in iel_to_orig_reg.items()}
    return elm_to_fr

def set_source_limits(cfg):
    geom = cfg.geometry
    br = geom.borehole.radius

    cfg_trans = cfg.transport_fullscale
    cfg_source = cfg_trans.source_params
    x_pos = cfg_source.source_ipos * (cfg_source.source_length + cfg_source.source_space)
    source_params = dict(
        source_y0=-2 * br,
        source_y1=2 * br,
        source_x0=x_pos,
        source_x1=x_pos + cfg_source.source_length,
    )
    return source_params


def compute_fields(cfg:dotdict, mesh:Mesh,  fr_map: Dict[int, Fracture]):
    """
    :param params: transport parameters dictionary
    :param mesh: GmshIO of the computational mesh (read only)
    :param fr_map: map ele id to the fracture (only for fracture 2d elements
    :return: el_ids:List[int], cond:List[float], cross:List[float]
    """
    cfg_geom = cfg.geometry
    cfg_trans = cfg.transport_fullscale

    #
    cfg_bulk_fields = cfg_trans.bulk_field_params
    conductivity = np.full( (len(mesh.elements),), float(cfg_bulk_fields.cond_min))
    cross_section = np.full( (len(mesh.elements),), float(1.0))
    porosity = np.full((len(mesh.elements),), 1.0)
    # Bulk fields
    el_slice_3d = mesh.el_dim_slice(3)
    bulk_cond, bulk_por = apply_fields.bulk_fields_mockup(cfg_geom, cfg_bulk_fields, mesh.el_barycenters()[el_slice_3d])
    conductivity[el_slice_3d] = bulk_cond
    porosity[el_slice_3d] = bulk_por
    logging.info(f"3D slice: {el_slice_3d}")
    c_min, c_max = np.min(conductivity), np.max(conductivity)
    logging.info(f"cond range: {c_min}, {c_max}")
    plots.plot_field(mesh.el_barycenters()[el_slice_3d], bulk_cond, file="conductivity_yz.pdf")
    plots.plot_field(mesh.el_barycenters()[el_slice_3d], bulk_por, file="porosity_yz.pdf")

    # Fracture
    cfg_fr_fields = cfg_trans.fr_field_params
    el_slice_2d = mesh.el_dim_slice(2)
    logging.info(f"2D slice: {el_slice_2d}")
    fr_cond, fr_cross, fr_por = apply_fields.fr_fields(cfg_fr_fields, mesh.elements[el_slice_2d], el_slice_2d.start, fr_map)
    conductivity[el_slice_2d] = fr_cond
    cross_section[el_slice_2d] = fr_cross
    porosity[el_slice_2d] = fr_por
    fields = dict(
        conductivity=conductivity,
        cross_section=cross_section,
        porosity=porosity
    )
    cond_file = mesh.write_fields("input_fields.msh2", fields)
    return cond_file


@report
@memoize
def fullscale_transport_mesh(cfg, cfg_mesh, seed):
    main_box_dimensions = cfg.geometry.box_dimensions
    fractures = mesh_tools.generate_fractures(cfg.fractures, main_box_dimensions, seed)
    return one_borehole(cfg.geometry, fractures, cfg_mesh), fractures