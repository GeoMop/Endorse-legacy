from typing import *
import os
import math
from bgem.gmsh import gmsh, field, options
import numpy as np
from .common import dotdict
import attrs


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

    borehole_distance = geom_dict.borehole_distance
    for i_shift in range(geom_dict.n_explicit_boreholes):
        laterals.append(lateral_tunnel_1.copy().translate([0, borehole_distance * i_shift, 0]))
        laterals.append(lateral_tunnel_1.copy().translate([0, -borehole_distance * i_shift, 0]))


    return main_tunnel_cylinder.fuse(*laterals)

def make_access_tunnels_simple(factory, geom_dict):
    r = geom_dict.main_tunnel.radius
    w = geom_dict.main_tunnel.width
    l = geom_dict.main_tunnel.length

    equivalent_height = w + r
    tunnel = factory.rectangle([equivalent_height, l]).rotate([0,1,0], math.pi/2).translate(
        [0, 0, equivalent_height / 2])
    return tunnel

def boreholes_full(factory, geom_dict):
    lateral_length = geom_dict.lateral_tunnel.length
    borehole_radius = geom_dict.borehole_radius
    borehole_length = geom_dict.borehole_length
    borehole_distance = geom_dict.borehole_distance

    b_1 = factory.cylinder(borehole_radius, axis=[borehole_length, 0, 0]).translate(
        [- borehole_length / 2, 0, 1.5 * borehole_radius])
    boreholes = [b_1]
    for i_shift in range(geom_dict.n_explicit_boreholes):
        boreholes.append(b_1.copy().translate([0, borehole_distance * i_shift, 0]))
        boreholes.append(b_1.copy().translate([0, -borehole_distance * i_shift, 0]))

    return factory.group(*boreholes)

def boreholes_simple(factory, geom_dict):
    lateral_length = geom_dict.lateral_tunnel.length
    borehole_radius = geom_dict.borehole_radius
    borehole_length = geom_dict.borehole_length
    borehole_distance = geom_dict.borehole_distance

    b_1 = factory.cylinder(borehole_radius, axis=[borehole_length, 0, 0]).translate(
        [-borehole_length / 2, 0, 1.5 * borehole_radius])

    # replace parallel boreholes by rectangles with matrix interface equivalent to the surface of 1.5 * borehole cylinder (approx boundary of EDZ)
    # interface accounts for both sides of the rectangle
    equivalent_height = borehole_radius * 1.5 * math.pi
    b_0 = factory.rectangle([borehole_length, equivalent_height]).rotate([1,0,0], math.pi/2).translate(
        [0, 0, equivalent_height / 2])
    b_2 = b_0.copy().translate([0, borehole_distance, 0])
    b_3 = b_0.translate([0, -borehole_distance, 0])


    return b_1, factory.group(b_2, b_3)

def basic_shapes(factory, geom_dict):
    box, sides = box_with_sides(factory, geom_dict.box_dimensions)
    borehole_length = geom_dict.borehole_length
    access_tunnels = make_access_tunnels(factory, geom_dict).translate([-borehole_length / 2, 0, 0])
    boreholes = boreholes_full(factory, geom_dict)
    tunnels = boreholes.copy().fuse(access_tunnels.copy())
    box_drilled = box.copy().cut(tunnels).set_region("box")
    return box_drilled, box, access_tunnels, boreholes

def basic_shapes_simple(factory, geom_dict):
    box, sides = box_with_sides(factory, geom_dict.box_dimensions)
    borehole_length = geom_dict.borehole_length
    access_tunnel_2d = make_access_tunnels_simple(factory, geom_dict).translate([-borehole_length / 2, 0, 0])
    borehole, boreholes_2d = boreholes_simple(factory, geom_dict)
    tunnels_2d = boreholes_2d.fuse(access_tunnel_2d.copy().cut(borehole.copy()))
    #tunnels = boreholes.copy().fragment(access_tunnels)

    box_fr, tunnels_2d_fr = factory.fragment(box, tunnels_2d)
    box_t2d = factory.group(box_fr, tunnels_2d_fr)
    #box_frag = gmsh.ObjectSet.group(box, tunnels_2d).copy().fragment(tunnels_2d.copy())
    #box_reg = gmsh.Region.get("box", 3)
    box_drilled = box_t2d.cut(borehole.copy())
    box_drilled.select_by_intersect(box).set_region("box")
    box_drilled.select_by_intersect(tunnels_2d).set_region("EDZ_2D")
    #box_drilled = factory.group(box_drilled, tunnels_2d)
    return box_drilled, box, borehole


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

def geom_meshing(factory, objects, mesh_file):
    factory.write_brep()
    factory.mesh_options.CharacteristicLengthMin = 0.1
    factory.mesh_options.CharacteristicLengthMax = 50
    factory.mesh_options.MinimumCirclePoints = 6
    factory.mesh_options.MinimumCurvePoints = 6
    #factory.mesh_options.Algorithm = options.Algorithm3d.MMG3D

    # mesh.Algorithm = options.Algorithm2d.MeshAdapt # produce some degenerated 2d elements on fracture boundaries ??
    # mesh.Algorithm = options.Algorithm2d.Delaunay
    # mesh.Algorithm = options.Algorithm2d.FrontalDelaunay

    factory.mesh_options.Algorithm = options.Algorithm3d.Delaunay
    #mesh.ToleranceInitialDelaunay = 0.01
    # mesh.ToleranceEdgeLength = fracture_mesh_step / 5
    #mesh.CharacteristicLengthFromPoints = True
    factory.mesh_options.CharacteristicLengthFromCurvature = False
    factory.mesh_options.CharacteristicLengthExtendFromBoundary = 2  # co se stane if 1
    #mesh.CharacteristicLengthMin = min_el_size
    #mesh.CharacteristicLengthMax = max_el_size

    #factory.keep_only(*objects)
    #factory.remove_duplicate_entities()
    factory.make_mesh(objects, dim=3)
    #factory.write_mesh(me gmsh.MeshFormat.msh2) # unfortunately GMSH only write in version 2 format for the extension 'msh2'
    factory.write_mesh(format=gmsh.MeshFormat.msh2)
    os.rename(factory.model_name + ".msh2", mesh_file)


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

    """
    cyl_dist = field.distance(center_line, sampling = 200) - 2
    ff = field.threshold(cyl_dist, lower_bound=(1, 0.3),  upper_bound=(50, 50))
    """
    """
    ff = field.attractor_aniso_curve(center_line,
                                     dist_range=(3, 10),
                                     h_normal=(0.3, 10),
                                     h_tangent=(3, 10),
                                     sampling=200)
    """

    #################################################################################

    center_line = factory.line([0,0,0], [cfg.borehole_length, 0, 0]).translate([-cfg.borehole_length / 2, 0, 1.5 * cfg.borehole_radius])

    bx, by, bz = cfg.box_dimensions
    n_sampling = int(cfg.borehole_length / 2)
    dist = field.distance(center_line, sampling = n_sampling)
    inner = field.geometric(dist, a=(cfg.borehole_radius, cfg.edz_mesh_step * 0.9), b=(cfg.edz_radius, cfg.edz_mesh_step))
    outer = field.polynomial(dist, a=(cfg.edz_radius, cfg.edz_mesh_step), b=(by/2, cfg.boundary_mesh_step), q=1.7)
    ff = field.maximum(inner, outer)

    factory.set_mesh_step_field(ff)
    geom_meshing(factory, [bulk], mesh_file)
    # factory.show()
    del factory
