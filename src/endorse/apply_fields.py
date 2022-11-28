import numpy as np

from .common import File
from .mesh_class import Mesh


def conductivity_mockup(cfg_geom, cfg_fields, output_mesh:Mesh):
    X, Y, Z = output_mesh.el_barycenters().T
    cond_file = "fine_conductivity.msh2"
    cond_max = float(cfg_fields.cond_max)
    cond_min = float(cfg_fields.cond_min)

    edz_r = cfg_geom.edz_radius # 2.5
    in_r = cfg_fields.inner_radius
    Z = Z - cfg_geom.borehole.z_pos
    # axis wit respect to EDZ radius
    Y_rel = Y / cfg_fields.h_axis
    Z_rel = Z / cfg_fields.v_axis

    # distance from center, 1== edz_radius
    distance = np.sqrt((Y_rel * Y_rel + Z_rel * Z_rel)) / (edz_r)

    theta = (1 - distance)/(1 - in_r)
    cond_field = np.minimum(cond_max, np.maximum(cond_min, np.exp(theta * np.log(cond_max) + (1-theta) * np.log(cond_min))))
    abs_dist = np.sqrt(Y * Y + Z * Z)
    cond_field[abs_dist < cfg_geom.borehole.radius] = 1e-18
    #print({(i+1):cond for i,cond in enumerate(cond_field)})
    output_mesh.write_fields(cond_file,
                            dict(conductivity=cond_field))
    return File(cond_file)


def bulk_fields_mockup(cfg_geom, cfg_fields, XYZ):
    X, Y, Z = XYZ.T

    edz_r = cfg_geom.edz_radius # 2.5
    in_r = cfg_fields.inner_radius
    Z = Z - cfg_geom.borehole.z_pos
    # axis wit respect to EDZ radius
    Y_rel = Y / cfg_fields.h_axis
    Z_rel = Z / cfg_fields.v_axis

    # distance from center, 1== edz_radius
    distance = np.sqrt((Y_rel * Y_rel + Z_rel * Z_rel)) / (edz_r)

    theta = (1 - distance)/(1 - in_r)
    cond_max = float(cfg_fields.cond_max)
    cond_min = float(cfg_fields.cond_min)
    cond_field = np.minimum(cond_max, np.maximum(cond_min, np.exp(theta * np.log(cond_max) + (1-theta) * np.log(cond_min))))
    #abs_dist = np.sqrt(Y * Y + Z * Z)
    #cond_field[abs_dist < cfg_geom.borehole.radius] = 1e-18

    por_max = float(cfg_fields.por_max)
    por_min = float(cfg_fields.por_min)
    por_field = np.minimum(por_max, np.maximum(por_min, np.exp(theta * np.log(por_max) + (1-theta) * np.log(por_min))))
    #cond_field[abs_dist < cfg_geom.borehole.radius] = 1e-18
    #print({(i+1):cond for i,cond in enumerate(cond_field)})

    return cond_field, por_field

viscosity = 1e-3
gravity_accel = 10
density = 1000
permeability_to_conductivity = viscosity * gravity_accel * density
def fr_fields(cfg_fr_fields, fr_elements, i_begin,  fr_map):
    """
    :param cfg_fr_fields:
    :param fr_elements: list of all 2d elements including the boundary elements
    :param i_begin:
    :param fr_map:
    :return:
    """
    fr_r = np.zeros(len(fr_elements))
    for iel, el in enumerate(fr_elements):
        try:
            fr = fr_map[iel + i_begin]
            fr_r[iel] = fr.r
        except KeyError:
            pass
    cross = float(cfg_fr_fields.apperture_per_size) * fr_r
    cond = permeability_to_conductivity/12 * cross * cross
    porosity = np.full_like(cond, 1.0)
    return cond, cross, porosity