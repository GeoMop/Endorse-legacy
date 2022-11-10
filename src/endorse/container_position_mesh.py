from typing import *
import os
import math
from bgem.gmsh import gmsh, field, options
import numpy as np
from .common import dotdict
from . import mesh_tools


def make_mesh(cfg:dotdict, fractures:List['Fracture'], i_pos:int, mesh_file: str):
    """
    Make mesh of given container position `i_pos`.
    BREP and mesh writen to given mesh_file derived files.
    The EDZ transport coordinate system is preserved.
        X - in direction of storage boreholes
        Y - perpendicular horizontal
        Z - vertical
        origin: center of the center borehole on the interface with lateral tunnel
    """

    # Radius of the homogenization kernel, approximately macro mesh step
    macro_mesh_step = 3
    b_cfg = cfg.borehole
    base, ext = os.path.splitext(os.path.basename(mesh_file))
    factory = gmsh.GeometryOCC(base, verbose=True)
    container_period = mesh_tools.container_period(cfg)
    box_shift = container_period/2 + mesh_tools.container_x_pos(cfg, i_pos)

    # TODO: homogenization x_size could be unrelated to the container size
    x_size = container_period + 2 * macro_mesh_step
    yz_size = 2 * (cfg.edz_radius + macro_mesh_step)
    box = factory.box([x_size, yz_size, yz_size]).translate([box_shift, 0, b_cfg.z_pos])
    bh = factory.cylinder(b_cfg.radius, axis=[b_cfg.length, 0, 0]).translate([0, 0, b_cfg.z_pos])
    domain = box.copy().fragment(bh.copy())
    outer = domain.select_by_intersect(box).set_region("outer")
    borehole = domain.select_by_intersect(bh).set_region("borehole")


    # TODO: mesh EDZ cylinder as well and implement region selection after meshing
    #edz = factory.cylinder(cfg.borehole_radius, axis=[cfg.borehole_length, 0, 0])

    factory.get_logger().start()
    #factory = gmsh.GeometryOCC(mesh_name)
    #gopt = options.Geometry()
    #gopt.Tolerance = 0.0001
    #gopt.ToleranceBoolean = 0.001
    factory.set_mesh_step_field(mesh_tools.edz_refinement_field(cfg, factory))
    factory.get_logger().stop()
    mesh_tools.edz_meshing(cfg, factory, [outer, borehole], mesh_file)
    # factory.show()
    del factory