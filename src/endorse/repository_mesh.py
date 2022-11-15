from typing import *
import os
import math
from bgem.gmsh import gmsh, options
import numpy as np
from endorse.common.common import dotdict
from . import mesh_tools

def create_fractures_rectangles(gmsh_geom, fractures, base_shape: 'ObjectSet'):
    # From given fracture date list 'fractures'.
    # transform the base_shape to fracture objects
    # fragment fractures by their intersections
    # return dict: fracture.region -> GMSHobject with corresponding fracture fragments
    if len(fractures) == 0:
        return []


    shapes = []
    for i, fr in enumerate(fractures):
        shape = base_shape.copy()
        print("fr: ", i, "tag: ", shape.dim_tags)
        shape = shape.scale([fr.rx, fr.ry, 1]) \
            .rotate(axis=fr.rotation_axis, angle=fr.rotation_angle) \
            .translate(fr.center).set_region(fr.region)
        shapes.append(shape)

    fracture_fragments = gmsh_geom.fragment(*shapes)
    return fracture_fragments


# @attrs.define
# class ThreeTunnelGeom:
#     main_tunnel: TunnelParams
#     lateral_tunnel: TunnelParams

def tunnel(factory, tunnel_dict):
    """
    A box with rounded "roof", basic box dimensions:
    hight= radius, width, length

    The crosscut is formed by the height x width rectangle with the dick segment roof.
    Result is translated to have [0,0,0] at the boundary of the floor rectangle.
    At the center of the 'width' side.
    """
    radius = tunnel_dict.radius
    height = tunnel_dict.height
    width = tunnel_dict.width
    length = tunnel_dict.length
    box = factory.box([width, length, height])
    z_shift = math.sqrt(radius * radius - 0.25 * width * width) - height / 2
    cylinder = factory.cylinder(radius, axis=[0, length, 0]).translate([0, -length / 2, -z_shift])
    roof = cylinder.intersect(box.copy().translate([0, 0, height]))
    return box.fuse(roof).translate([0,0,+height / 2])

def box_with_sides(factory, dimensions):
    """
    Make a box and dictionary of its sides named: 'side_[xyz][01]'
    :return: box, sides_dict
    """
    box = factory.box(dimensions).set_region("box")
    side_z = factory.rectangle([dimensions[0], dimensions[1]])
    side_y = factory.rectangle([dimensions[0], dimensions[2]])
    side_x = factory.rectangle([dimensions[2], dimensions[1]])
    sides = dict(
        side_z0=side_z.copy().translate([0, 0, -dimensions[2] / 2]),
        side_z1=side_z.copy().translate([0, 0, +dimensions[2] / 2]),
        side_y0=side_y.copy().translate([0, 0, -dimensions[1] / 2]).rotate([-1, 0, 0], np.pi / 2),
        side_y1=side_y.copy().translate([0, 0, +dimensions[1] / 2]).rotate([-1, 0, 0], np.pi / 2),
        side_x0=side_x.copy().translate([0, 0, -dimensions[0] / 2]).rotate([0, 1, 0], np.pi / 2),
        side_x1=side_x.copy().translate([0, 0, +dimensions[0] / 2]).rotate([0, 1, 0], np.pi / 2)
    )
    for name, side in sides.items():
        side.modify_regions(name)
    return box, sides

def make_access_tunnels(factory, geom_dict):

    lateral_length = geom_dict.lateral_tunnel.length
    main_tunnel_cylinder = tunnel(
        factory, geom_dict.main_tunnel,
        ).translate([-lateral_length, 0, 0])
    lateral_tunnel_1 = tunnel(
        factory, geom_dict.lateral_tunnel,
    ).rotate([0,0,1], math.pi/2).translate([-lateral_length/2, 0, 0])

    laterals = [lateral_tunnel_1]

    borehole_distance = geom_dict.borehole.y_spacing
    for i_shift in range(geom_dict.borehole.n_explicit):
        laterals.append(lateral_tunnel_1.copy().translate([0, borehole_distance * i_shift, 0]))
        laterals.append(lateral_tunnel_1.copy().translate([0, -borehole_distance * i_shift, 0]))


    return main_tunnel_cylinder.fuse(*laterals)

def boreholes_full(factory, geom_dict):
    lateral_length = geom_dict.lateral_tunnel.length
    b_cfg = geom_dict.borehole
    borehole_radius = b_cfg.radius
    borehole_length = b_cfg.length
    borehole_distance = b_cfg.y_spacing

    b_1 = factory.cylinder(borehole_radius, axis=[borehole_length, 0, 0])
    boreholes = [b_1]
    for i_shift in range(geom_dict.borehole.n_explicit):
        boreholes.append(b_1.copy().translate([0, borehole_distance * i_shift, 0]))
        boreholes.append(b_1.copy().translate([0, -borehole_distance * i_shift, 0]))

    return factory.group(*boreholes)

def basic_shapes(factory, geom_dict):
    bh_length = geom_dict.borehole.length
    bh_z_pos = geom_dict.borehole.z_pos

    box, sides = box_with_sides(factory, geom_dict.box_dimensions)
    box = box.translate([bh_length / 2, 0, 0])
    access_tunnels = make_access_tunnels(factory, geom_dict) #.translate([-bh_length / 2, 0, 0])
    boreholes = boreholes_full(factory, geom_dict).translate([0, 0, bh_z_pos])
    tunnels = boreholes.copy().fuse(access_tunnels.copy())
    box_drilled = box.copy().cut(tunnels).set_region("box")
    return box_drilled, box, access_tunnels, boreholes


def make_geometry(factory, geom_dict, fractures):
    box_drilled, box, access_tunnels, boreholes = basic_shapes(factory, geom_dict)
    #box_drilled, box, tunnels = basic_shapes_simple(factory, geom_dict)

    fractures = create_fractures_rectangles(factory, fractures, factory.rectangle())
    fractures_group = factory.group(*fractures).intersect(box_drilled)

    #b_rec = box_drilled.get_boundary()#.set_region(".sides")

    box_fr, fractures_fr = factory.fragment(box_drilled, fractures_group)
    fractures_fr.mesh_step(geom_dict.fracture_mesh_step) #.set_region("fractures")

    b_box_fr = box_fr.get_boundary().split_by_dimension()[2]
    b_fractures_fr = fractures_fr.get_boundary().split_by_dimension()[1]

    # select outer boundary
    boundary_mesh_step = geom_dict.boundary_mesh_step
    b_box = b_box_fr.select_by_intersect(box.get_boundary().copy()).set_region(".box_outer").mesh_step(boundary_mesh_step)
    b_fractures = b_fractures_fr.select_by_intersect(box.get_boundary().copy()).set_region(".fr_outer").mesh_step(boundary_mesh_step)

    # select inner boreholes boundary
    boreholes_step = geom_dict.boreholes_mesh_step
    select = boreholes.get_boundary().copy()
    b_box_boreholes = b_box_fr.select_by_intersect(select)\
                  .set_region(".box_boreholes").mesh_step(boreholes_step)
    b_fr_boreholes = b_fractures_fr.select_by_intersect(select)\
                 .set_region(".fr_boreholes").mesh_step(boreholes_step)

    tunnel_mesh_step = geom_dict.main_tunnel_mesh_step
    select = access_tunnels.get_boundary().copy()
    b_box_tunnel = b_box_fr.select_by_intersect(select)\
                  .set_region(".box_tunnel").mesh_step(tunnel_mesh_step)
    b_fr_tunnel = b_fractures_fr.select_by_intersect(select)\
                  .set_region(".fr_tunnel").mesh_step(tunnel_mesh_step)


    boundary = factory.group(b_box, b_fractures,
                             b_box_boreholes, b_fr_boreholes,
                             b_box_tunnel, b_fr_tunnel)
    bulk_geom = factory.group(box_fr, fractures_fr, boundary)
    edz_refined = factory.group(b_box_boreholes, b_fr_boreholes, b_box_tunnel, b_fr_tunnel)
    #boundary = factory.group(b_box)

    # Following makes some mesing issues:
    #factory.group(b_box_inner, b_fr_inner).mesh_step(geom_dict['main_tunnel_mesh_step'])
    #boundary.select_by_intersect(boreholes.get_boundary()).mesh_step(geom_dict['boreholes_mesh_step'])

    return bulk_geom, edz_refined


def make_mesh(cfg:dotdict, fractures:List['Fracture'], mesh_file: str):
    """
    :param cfg: repository mesh configuration cfg.repository_mesh
    :param fractures:  generated fractures
    :param mesh_file:
    :return:
    """
    base, ext = os.path.splitext(os.path.basename(mesh_file))
    factory = gmsh.GeometryOCC(base, verbose=True)
    factory.get_logger().start()
    #factory = gmsh.GeometryOCC(mesh_name)
    gopt = options.Geometry()
    gopt.Tolerance = 0.0001
    gopt.ToleranceBoolean = 0.001

    bulk, refined = make_geometry(factory, cfg, fractures)


    factory.set_mesh_step_field(mesh_tools.edz_refinement_field(cfg, factory))
    mesh_tools.edz_meshing(cfg, factory, [bulk], mesh_file)
    # factory.show()
    del factory
