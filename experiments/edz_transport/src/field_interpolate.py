from typing import *
from bgem.gmsh.gmsh_io import GmshIO
from bgem.stochastic.fracture import Fracture
import numpy as np

def compute_fields(params, mesh,  fr_map: Dict[int, Fracture]):
    """
    :param params: transport parameters dictionary
    :param mesh: GmshIO of the computational mesh (read only)
    :param fr_map: map ele id to the fracture (only for fracture 2d elements
    :return: el_ids:List[int], cond:List[float], cross:List[float]
    """
    el_ids = []
    cond = []
    cross = []
    bulk_cond = float(params['bulk_conductivity'])
    apperture_per_size = float(params['apperture_per_size'])
    for id, el in mesh.elements.items():
        (type, tags, nodes) = el
        i_reg, i_shape = tags
        conductivity = None
        if type == 4: # 3d
            # TODO: interpolate from 2D VTK results
            conductivity = bulk_cond
            cross_section = 1.0
        elif type == 2: # 2d
            if  id in fr_map: # fracture element
                fr = fr_map[id]
                cross_section = apperture_per_size * fr.r
                conductivity = 1/12 * cross_section * cross_section * 1000 * 10 * 1000
            elif mesh.physical["EDZ_2D"] == (i_reg, 2): # simplified edz element
                # effective thickness of EDZ about 1m
                cross_section = 1
                # effective conductivity (TODO: averaging)
                # One order of magnitude above bulk.
                conductivity = bulk_cond * 10

        if conductivity is not None:
            el_ids.append(id)
            cond.append(conductivity)
            cross.append(cross_section)
    return el_ids, cond, cross

def field_interpolate(cfg, el_fr_map):
    input_mesh = "random_fractures_healed.msh"
    mesh = GmshIO(filename=input_mesh)
    ele_ids, cond, cross = compute_fields(cfg['transport_params'], mesh, el_fr_map)
    fields = {
        "conductivity" : np.array(cond)[:, None],
        "cross_section" : np.array(cross)[:, None]
    }
    mesh.write_ascii("input_fields.msh")
    mesh.write_fields("input_fields.msh", ele_ids, fields)

"""
Different number of entities: 157203, computation needs 157407.
.. missing values for fractures !!

mesh:
bulk tetrahedra: 155774
EDZ: 204
fractures: 1429 
"""