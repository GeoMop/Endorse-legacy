import logging
import os
from typing import *

import numpy as np
from endorse.mesh.repository_mesh import fullscale_transport_mesh_3d, fullscale_transport_mesh_2d

from . import common
from .common import dotdict, File, report, memoize
from .mesh_class import Mesh
from . import apply_fields
from . import plots
from . import flow123d_inputs_path
from .indicator import indicators, IndicatorFn
from bgem.stochastic.fracture import Fracture, Population


# def input_files(cfg):
#     return [
#         cfg_tr_full.piezo_head_input_file,
#         cfg_tr_full.conc_flux_file
#     ]

#
# def find_nearest(array, value):
#     array = np.asarray(array)
#     idx = (np.abs(array - value)).argmin()
#     return idx, array[idx]

def fullscale_transport(cfg_path, seed):
    cfg = common.load_config(cfg_path)
    return transport_run(cfg, seed)

def transport_2d(cfg, seed):
    """
    1. apply conouctivity to given mesh:
       - on borehole neighbourhood, select elements
       - calculate barycenters
       - apply conductivity
       - write the field
    2. substitute source term space distribution
    3. return necessary files
    """
    #files = input_files(cfg.transport_fullscale)
    cfg_basedir = cfg._config_root_dir
    cfg_fine = cfg.transport_fullscale
    large_model = File(os.path.join(cfg_basedir, cfg_fine.piezo_head_input_file))
    #conc_flux = File(os.path.join(cfg_basedir, cfg_fine.conc_flux_file))
    #plots.plot_source(conc_flux)

    box = cfg.geometry.box_dimensions
    box = [box[0], box[2], 0]
    fr_pop = Population.initialize_2d( cfg.fractures.population, box)

    full_mesh_file, fractures, n_large = fullscale_transport_mesh_2d(cfg_fine, fr_pop, seed)

    full_mesh = Mesh.load_mesh(full_mesh_file, heal_tol=1e-4)
    el_to_ifr = fracture_map(full_mesh, fractures, n_large)
    # mesh_modified_file = full_mesh.write_fields("mesh_modified.msh2")
    # mesh_modified = Mesh.load_mesh(mesh_modified_file)

    input_fields_file, est_velocity = compute_fields(cfg, full_mesh, el_to_ifr, fractures)

    # input_fields_file = compute_fields(cfg, full_mesh, el_to_fr)
    params = cfg_fine.copy()

    # estimate times
    #bulk_vel_est, fr_vel_est = est_velocity
    #end_time = (50 / bulk_vel_est + 50 / fr_vel_est)
    #dt = 0.5 / bulk_vel_est
    # convert to years

    #end_time = end_time / common.year
    #dt = dt / common.year

    #end_time = 10 * dt
    new_params = dict(
        mesh_file=input_fields_file,
        piezo_head_input_file=large_model,
        #conc_flux_file=conc_flux,
        input_fields_file = input_fields_file,
        dg_penalty = cfg_fine.dg_penalty,
        end_time_years = cfg_fine.end_time,
        trans_solver__a_tol= cfg_fine.trans_solver__a_tol,
        trans_solver__r_tol= cfg_fine.trans_solver__r_tol

        #max_time_step = dt,
        #output_step = 10 * dt
    )
    params.update(new_params)
    params.update(set_source_limits(cfg))
    template = flow123d_inputs_path.joinpath(cfg_fine.input_template)
    fo = common.call_flow(cfg.flow_env, template, params)
    return get_indicator(cfg, fo)


def transport_run(cfg, seed):
    """
    1. apply conouctivity to given mesh:
       - on borehole neighbourhood, select elements
       - calculate barycenters
       - apply conductivity
       - write the field
    2. substitute source term space distribution
    3. return necessary files
    """
    #files = input_files(cfg.transport_fullscale)
    cfg_basedir = cfg._config_root_dir
    cfg_fine = cfg.transport_fullscale
    large_model = File(os.path.join(cfg_basedir, cfg_fine.piezo_head_input_file))
    #conc_flux = File(os.path.join(cfg_basedir, cfg_fine.conc_flux_file))
    #plots.plot_source(conc_flux)

    fr_pop = Population.initialize_3d( cfg.fractures.population, cfg.geometry.box_dimensions)

    full_mesh_file, fractures, n_large = fullscale_transport_mesh_3d(cfg_fine, fr_pop, seed)

    full_mesh = Mesh.load_mesh(full_mesh_file, heal_tol=1e-4)
    el_to_ifr = fracture_map(full_mesh, fractures, n_large)
    # mesh_modified_file = full_mesh.write_fields("mesh_modified.msh2")
    # mesh_modified = Mesh.load_mesh(mesh_modified_file)

    input_fields_file, est_velocity = compute_fields(cfg, full_mesh, el_to_ifr, fractures)

    # input_fields_file = compute_fields(cfg, full_mesh, el_to_fr)
    params = cfg_fine.copy()

    # estimate times
    #bulk_vel_est, fr_vel_est = est_velocity
    #end_time = (50 / bulk_vel_est + 50 / fr_vel_est)
    #dt = 0.5 / bulk_vel_est
    # convert to years

    #end_time = end_time / common.year
    #dt = dt / common.year

    #end_time = 10 * dt
    new_params = dict(
        mesh_file=input_fields_file,
        piezo_head_input_file=large_model,
        #conc_flux_file=conc_flux,
        input_fields_file = input_fields_file,
        dg_penalty = cfg_fine.dg_penalty,
        end_time_years = cfg_fine.end_time,
        trans_solver__a_tol= cfg_fine.trans_solver__a_tol,
        trans_solver__r_tol= cfg_fine.trans_solver__r_tol

        #max_time_step = dt,
        #output_step = 10 * dt
    )
    params.update(new_params)
    params.update(set_source_limits(cfg))
    template = flow123d_inputs_path.joinpath(cfg_fine.input_template)
    fo = common.call_flow(cfg.flow_env, template, params)
    return get_indicator(cfg, fo)

@report
def get_indicator(cfg, fo):
    cfg_fine = cfg.transport_fullscale
    z_dim = 0.9 * 0.5 * cfg.geometry.box_dimensions[2]
    z_shift = cfg.geometry.borehole.z_pos
    z_cuts = (z_shift - z_dim, z_shift + z_dim)
    inds = indicators(fo.solute.spatial_file, f"{cfg_fine.conc_name}_conc", z_cuts)
    plots.plot_indicators(inds)
    itime = IndicatorFn.common_max_time(inds)  # not splined version, need slice data
    #plots.plot_slices(fo.solute.spatial_file, f"{cfg_fine.conc_name}_conc", z_cuts, [itime-1, itime, itime+1])
    ind_time_max = [ind.time_max()[1] for ind in inds]
    return ind_time_max

@report
def fracture_map(mesh, fractures, n_large) -> Dict[int, Fracture]:
    """
    - join all fracture regions into single "fractures" region
    - return dictionary mapping element idx to fracture
    :param mesh:
    :param fractures:
    :return:
    """
    own_name_to_id = {fr.region.name: fr.region.id for fr in fractures}
    own_to_gmsh_id = {own_name_to_id[name]: gmsh_id for name,(gmsh_id, dim) in mesh.gmsh_io.physical.items() if name in own_name_to_id}

    max_reg = max( [gmsh_id for gmsh_id, dim in mesh.gmsh_io.physical.values()] )
    small_reg_id = max_reg + 1
    large_reg_id = max_reg + 2
    large_reg_map = {(own_to_gmsh_id[fr.region.id], 1): (large_reg_id, 1, "fr_large")  for fr in fractures[:n_large] if fr.region.id in own_to_gmsh_id}
    small_reg_map = {(own_to_gmsh_id[fr.region.id], 1): (small_reg_id, 1, "fr_small")  for fr in fractures[n_large:] if fr.region.id in own_to_gmsh_id}
    new_reg_map = large_reg_map
    new_reg_map.update(small_reg_map)

    # if do_heal:
    #     hm.heal_mesh(gamma_tol=0.01)
    #     hm.move_all(geom_dict["shift_vec"])
    #     elm_to_orig_reg = hm.map_regions(new_reg_map)
    #     hm.stats_to_yaml(mesh_name + "_heal_stats.yaml")
    #     assert hm.healed_mesh_name == mesh_healed
    #     hm.write()
    # else:
    iel_to_orig_reg = mesh.map_regions(new_reg_map)

    reg_to_ifr = {own_to_gmsh_id[fr.region.id]: ifr for ifr, fr in enumerate(fractures) if fr.region.id in own_to_gmsh_id }
    elm_to_ifr = {el_idx: reg_to_ifr[reg_id] for el_idx, (reg_id, dim) in iel_to_orig_reg.items()}
    return elm_to_ifr

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

@report
def compute_fields(cfg:dotdict, mesh:Mesh,  fr_map: Dict[int, int], fractures:List[Fracture]):
    """
    :param params: transport parameters dictionary
    :param mesh: GmshIO of the computational mesh (read only)
    :param fr_map: map ele id to the fracture (only for fracture 2d elements
    :return: el_ids:List[int], cond:List[float], cross:List[float]
    """
    cfg_geom = cfg.geometry
    cfg_trans = cfg.transport_fullscale


    cfg_bulk_fields = cfg_trans.bulk_field_params

    conductivity = np.full( (len(mesh.elements),), float(cfg_bulk_fields.cond_min))
    cross_section = np.full( (len(mesh.elements),), float(1.0))
    porosity = np.full((len(mesh.elements),), 1.0)
    # Bulk fields
    el_slice_2d = mesh.el_dim_slice(2)
    bulk_cond, bulk_por = apply_fields.bulk_fields_mockup(cfg_geom, cfg_bulk_fields, mesh.el_barycenters()[el_slice_2d])
    conductivity[el_slice_2d] = bulk_cond
    porosity[el_slice_2d] = bulk_por
    logging.info(f"2D slice: {el_slice_2d}")
    c_min, c_max = np.min(conductivity), np.max(conductivity)
    logging.info(f"cond range: {c_min}, {c_max}")
    plots.plot_field(mesh.el_barycenters()[el_slice_2d], bulk_cond, file="conductivity_yz.pdf")
    plots.plot_field(mesh.el_barycenters()[el_slice_2d], bulk_por, file="porosity_yz.pdf")

    # Fracture
    cfg_fr = cfg_trans.fractures
    cfg_fr_fields = cfg_trans.fr_field_params
    el_slice_1d = mesh.el_dim_slice(1)
    logging.info(f"2D slice: {el_slice_1d}")
    i_default = len(fractures)
    fr_map_slice = [fr_map.get(i, i_default) for i in range(el_slice_1d.start, el_slice_1d.stop)]
    fr_cond, fr_cross, fr_por = apply_fields.fr_fields_repo(cfg_fr, cfg_fr_fields,
                                                            mesh.elements[el_slice_1d], fr_map_slice, fractures)
    conductivity[el_slice_1d] = fr_cond
    cross_section[el_slice_1d] = fr_cross
    porosity[el_slice_1d] = fr_por
    fields = dict(
        conductivity=conductivity,
        cross_section=cross_section,
        porosity=porosity
    )
    cond_file = mesh.write_fields("input_fields.msh2", fields)

    # estimate velocities on bulk and fracture
    # for cond range 1e-13 - 1e-9 and porosity about 1, we have velocity 1e-16  to 5.5e-10
    # i.e velocity about the order of conductivity or one order less
    # for fracture, cond range:

    pos_fr = fr_cond > 0
    est_velocity = (np.quantile(bulk_cond, 0.4)/10, np.quantile(fr_cond[pos_fr],  0.4))
    return cond_file, est_velocity

# def transport_observe_points(cfg):
#     cfg_geom = cfg.geometry
#     lx, ly, lz = cfg_geom.box_dimensions
#     np.
#     with "observe_points.csv":
#
# observe_points:
# - [0, 0.1, 0]
# - {point: [0.55, 0.55, 0], snap_region: 1d_lower}
# - {point: [0.7, 0.8, 0], snap_region: 1d_upper}
