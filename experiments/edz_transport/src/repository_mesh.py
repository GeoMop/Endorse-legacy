from bgem.gmsh import gmsh, field, options
import numpy as np
import math
import os
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

def tunnel(factory, radius, width, length):
    """
    A box with rounded "roof", basic box dimensions:
    hight= radius, width, length

    The crosscut is formed by the height x width rectangle with the dick segment roof.
    Result is translated to have [0,0,0] at the boundary of the floor rectangle.
    At the center of the 'width' side.
    """
    height = radius
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

    lateral_length = geom_dict["lateral_tunnel_length"]
    main_tunnel_cylinder = tunnel(
        factory,
        geom_dict["main_tunnel_radius"],
        geom_dict["main_tunnel_width"],
        geom_dict["main_tunnel_length"]
        ).translate([-lateral_length, 0, 0])
    lateral_tunnel_1 = tunnel(
        factory,
        geom_dict["lateral_tunnel_radius"],
        geom_dict["lateral_tunnel_width"],
        lateral_length
    ).rotate([0,0,1], math.pi/2).translate([-lateral_length/2, 0, 0])

    borehole_distance = geom_dict["borehole_distance"]
    lateral_tunnel_2 = lateral_tunnel_1.copy().translate([0, borehole_distance, 0])
    lateral_tunnel_3 = lateral_tunnel_1.copy().translate([0, -borehole_distance, 0])


    return main_tunnel_cylinder.fuse(lateral_tunnel_1, lateral_tunnel_2, lateral_tunnel_3)

def make_access_tunnels_simple(factory, geom_dict):
    r = geom_dict["main_tunnel_radius"]
    w = geom_dict["main_tunnel_width"]
    l = geom_dict["main_tunnel_length"]

    equivalent_height = w + r
    tunnel = factory.rectangle([equivalent_height, l]).rotate([0,1,0], math.pi/2).translate(
        [0, 0, equivalent_height / 2])
    return tunnel

def boreholes_full(factory, geom_dict):
    lateral_length = geom_dict["lateral_tunnel_length"]
    borehole_radius = geom_dict["borehole_radius"]
    borehole_length = geom_dict["borehole_length"]
    borehole_distance = geom_dict["borehole_distance"]

    b_1 = factory.cylinder(borehole_radius, axis=[borehole_length, 0, 0]).translate(
        [- borehole_length / 2, 0, 1.5 * borehole_radius])
    b_2 = b_1.copy().translate([0, borehole_distance, 0])
    b_3 = b_1.copy().translate([0, -borehole_distance, 0])
    return factory.group(b_1, b_2, b_3)

def boreholes_simple(factory, geom_dict):
    lateral_length = geom_dict["lateral_tunnel_length"]
    borehole_radius = geom_dict["borehole_radius"]
    borehole_length = geom_dict["borehole_length"]
    borehole_distance = geom_dict["borehole_distance"]

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
    box, sides = box_with_sides(factory, geom_dict["box_dimensions"])
    borehole_length = geom_dict["borehole_length"]
    access_tunnels = make_access_tunnels(factory, geom_dict).translate([-borehole_length / 2, 0, 0])
    boreholes = boreholes_full(factory, geom_dict)
    tunnels = boreholes.copy().fuse(access_tunnels)
    box_drilled = box.copy().cut(tunnels).set_region("box")
    return box_drilled, box, tunnels

def basic_shapes_simple(factory, geom_dict):
    box, sides = box_with_sides(factory, geom_dict["box_dimensions"])
    borehole_length = geom_dict["borehole_length"]
    access_tunnel_2d = make_access_tunnels_simple(factory, geom_dict).translate([-borehole_length / 2, 0, 0])
    borehole, boreholes_2d = boreholes_simple(factory, geom_dict)
    tunnels_2d = boreholes_2d.fuse(access_tunnel_2d.copy().cut(borehole.copy()))
    #tunnels = boreholes.copy().fragment(access_tunnels)

    box_drilled = box.copy().fragment(tunnels_2d).cut(borehole).set_region("box")
    tunnels_2d.set_region("EDZ_2D")
    box_drilled = factory.group(box_drilled, tunnels_2d)
    return box_drilled, box, borehole


def make_geometry(factory, geom_dict, fractures):
    #box_drilled, box, tunnels = basic_shapes(factory, geom_dict)
    box_drilled, box, tunnels = basic_shapes_simple(factory, geom_dict)

    fractures = create_fractures_rectangles(factory, fractures, factory.rectangle())
    fractures_group = factory.group(*fractures).intersect(box_drilled)

    #b_rec = box_drilled.get_boundary()#.set_region(".sides")

    box_fr, fractures_fr = factory.fragment(box_drilled, fractures_group)
    fractures_fr.mesh_step(geom_dict['fracture_mesh_step']) #.set_region("fractures")

    bulk_geom = factory.group(box_fr, fractures_fr)
    b_box_fr = box_fr.get_boundary().split_by_dimension()[2]
    b_fractures_fr = fractures_fr.get_boundary().split_by_dimension()[1]
    boundary_mesh_step = geom_dict['boundary_mesh_step']
    b_box = b_box_fr.select_by_intersect(box.get_boundary().copy()).set_region(".box_outer").mesh_step(boundary_mesh_step)
    b_fractures = b_fractures_fr.select_by_intersect(box.get_boundary().copy()).set_region(".fr_outer").mesh_step(boundary_mesh_step)
    b_box_inner = b_box_fr.select_by_intersect(tunnels.get_boundary().copy()).set_region(".box_inner").mesh_step(geom_dict['boreholes_mesh_step'])
    b_fr_inner = b_fractures_fr.select_by_intersect(tunnels.get_boundary().copy()).set_region(".fr_inner").mesh_step(geom_dict['boreholes_mesh_step'])
    boundary = factory.group(b_box, b_fractures, b_box_inner, b_fr_inner)
    #boundary = factory.group(b_box)

    # Following makes some mesing issues:
    #factory.group(b_box_inner, b_fr_inner).mesh_step(geom_dict['main_tunnel_mesh_step'])
    #boundary.select_by_intersect(boreholes.get_boundary()).mesh_step(geom_dict['boreholes_mesh_step'])

    return bulk_geom, boundary

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


def make_mesh(geom_dict, fractures, mesh_file):
    base, ext = os.path.splitext(os.path.basename(mesh_file))
    factory = gmsh.GeometryOCC(base, verbose=True)
    #factory = gmsh.GeometryOCC(mesh_name)
    gopt = options.Geometry()
    gopt.Tolerance = 0.0001
    gopt.ToleranceBoolean = 0.001

    bulk, boundary = make_geometry(factory, geom_dict, fractures)

    #a, b, c, d = 2, 3, 5, 15
    #f_borehole = field.threshold(field.distance_surfaces(b_borehole_fr.tags, 2 * bl / (2 * b)), lower_bound=(a, b), upper_bound=(c, d))
    #f_borehole = field.threshold(field.distance_surfaces(b_rec_fr.tags, 2 * bl / (2 * b)), lower_bound=(a, b), upper_bound=(c, d))
    borehole_radius = geom_dict["borehole_radius"]
    borehole_length = geom_dict["borehole_length"]
    center_line = factory.line([0,0,0], [borehole_length, 0, 0]).translate([-borehole_length / 2, 0, 1.5 * borehole_radius])
    """
    cyl_dist = field.distance(center_line, sampling = 200) - 2
    ff = field.threshold(cyl_dist, lower_bound=(1, 0.3),  upper_bound=(50, 50))
    """

    ff = field.attractor_aniso_curve(center_line,
                                     dist_range=(3, 10),
                                     h_normal=(0.3, 10),
                                     h_tangent=(3, 10),
                                     sampling=200)
    #################################################################################
    '''
    
   #block_fr = block.fragment(fractures_group)
    b_tmp = block.get_boundary()
    b_block_fr = b_rec_fr.select_by_intersect(b_tmp).set_region(".main_tunnel")

    a, b, c, d = 3, 5, 5, 15
    f_block = field.threshold(field.distance_surfaces(b_block_fr.tags, 2 * bl / (2 * b)), lower_bound=(a, b),
                        upper_bound=(c, d))

    #################################################################################

    ff = field.minimum(f_borehole, f_block)
    '''


    ff = field.constant(10)
    factory.set_mesh_step_field(ff)
    geom_meshing(factory, [bulk, boundary], mesh_file)
    # factory.show()
    del factory
