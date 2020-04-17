import numpy as np
import math
import os

from bgem.gmsh import gmsh
from bgem.gmsh import options
from bgem.gmsh import field
import process

def shift(radius, width):
    return math.sqrt(radius * radius - 0.25 * width * width) - radius / 2

def make_mesh(config_dict, fractures, mesh_name, mesh_file):
    geom = config_dict["geometry"]
    fracture_mesh_step = geom['fracture_mesh_step']
    dimensions = geom["box_dimensions"]
    # # well_z0, well_z1 = geom["well_openning"]
    # # well_length = well_z1 - well_z0
    # # well_r = geom["well_effective_radius"]
    # well_dist = geom["well_distance"]

    mtr = geom["main_tunnel_radius"]
    mtw = geom["main_tunnel_width"]
    mtl = geom["main_tunnel_length"]

    str = geom["lateral_tunnel_radius"]
    stw = geom["lateral_tunnel_width"]
    stl = geom["lateral_tunnel_length"]
    # small_tunnel_passage_width = geom["lateral_tunnel_passage_width"]

    br = geom["borehole_radius"]
    bl = geom["borehole_length"]
    bd = geom["borehole_distance"]

    factory = gmsh.GeometryOCC(mesh_name, verbose=True)
    gopt = options.Geometry()
    gopt.Tolerance = 0.0001
    gopt.ToleranceBoolean = 0.001
    # gopt.MatchMeshTolerance = 1e-1

    # Main box
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

    b_box = box.get_boundary().copy()

    # # two vertical cut-off wells, just permeable part
    # left_center = [-well_dist/2, 0, 0]
    # right_center = [+well_dist/2, 0, 0]

    # longitudinal part of main tunnel
    main_tunnel_block = factory.box([mtw, mtr, mtl])

    main_tunnel_block_tmp = main_tunnel_block.copy().translate([0, mtr, 0])
    main_tunnel_cylinder_tmp = factory.cylinder(mtr, axis=[0, 0, mtl]).translate([0, -shift(mtr, mtw), -mtl/2])
    main_tunnel_cylindrical_part = main_tunnel_block_tmp.intersect(main_tunnel_cylinder_tmp)

    # lateral part of main tunnel
    small_tunnel_block = factory.box([stl, str, stw])

    small_tunnel_block_tmp = small_tunnel_block.copy().translate([0, str, 0])
    small_tunnel_cylinder_tmp = factory.cylinder(str, axis=[stl, 0, 0]).translate([-stl/2, -shift(str, stw), 0])
    small_tunnel_cylindrical_part = small_tunnel_block_tmp.intersect(small_tunnel_cylinder_tmp)

    # lateral_tunnel_part_1 = small_tunnel_cylindrical_part.fuse(small_tunnel_block).translate([stl/2, 0, 0])
    # lateral_tunnel_part_2 = lateral_tunnel_part_1.copy().translate([0, 0, bd])
    # # lateral_tunnel_part = factory.group(lateral_tunnel_part_1, lateral_tunnel_part_2)

    # # main tunnel with lateral parts
    # main_tunnel = main_tunnel_block.fuse(main_tunnel_cylindrical_part, lateral_tunnel_part_1, lateral_tunnel_part_2)
    #     # .translate([-borehole_length/2, 0, 0])

    main_tunnel = main_tunnel_block.fuse(main_tunnel_cylindrical_part)

    # horizontal boreholes
    borehole1 = factory.cylinder(br, axis=[bl, 0, 0]).translate([stl, 0, 0])
    borehole2 = borehole1.copy().translate([0, 0, bd])
    borehole = borehole1.fuse(borehole2)

    # # right_well = factory.group(right_well_1, right_well_2).translate([-30, 0, 0])
    # boreholes = right_well_1.fuse(right_well_2).translate([-30, 0, 0])

    b_main_tunnel = main_tunnel.get_boundary()
    b_borehole = borehole.get_boundary()
    # b_borehole2 = borehole2.get_boundary()

    print("n fractures:", len(fractures))
    fractures = process.create_fractures_rectangles(factory, fractures, factory.rectangle())
    #fractures = create_fractures_polygons(factory, fractures)
    fractures_group = factory.group(*fractures)
    #fractures_group = fractures_group.remove_small_mass(fracture_mesh_step * fracture_mesh_step / 10)

    # drilled box and its boundary
    box_drilled = box.cut(main_tunnel, borehole)

    # fractures, fragmented, fractures boundary
    print("cut fractures by box without wells")
    fractures_group = fractures_group.intersect(box_drilled.copy())
    print("fragment fractures")
    box_fr, fractures_fr = factory.fragment(box_drilled, fractures_group)
    print("finish geometry")
    b_box_fr = box_fr.get_boundary()
    b_main_tunnel_fr = b_box_fr.select_by_intersect(b_main_tunnel).set_region(".main_tunnel")
    b_borehole_fr = b_box_fr.select_by_intersect(b_borehole).set_region(".borehole")
    # b_right_r2 = b_box_fr.select_by_intersect(b_borehole2).set_region(".borehole2")

    box_all = []
    for name, side_tool in sides.items():
        isec = b_box_fr.select_by_intersect(side_tool)
        box_all.append(isec.modify_regions("." + name))
    box_all.extend([box_fr, b_main_tunnel_fr, b_borehole_fr])

    b_fractures = factory.group(*fractures_fr.get_boundary_per_region())
    b_fractures_box = b_fractures.select_by_intersect(b_box).modify_regions("{}_box")
    b_fr_main_tunnel = b_fractures.select_by_intersect(b_main_tunnel).modify_regions("{}_main_tunnel")
    b_fr_borehole = b_fractures.select_by_intersect(b_borehole).modify_regions("{}_borehole")
    # b_fr_bor_2 = b_fractures.select_by_intersect(b_borehole2).modify_regions("{}_borehole2")

    b_fractures = factory.group(b_fr_main_tunnel, b_fr_borehole, b_fractures_box)
    mesh_groups = [*box_all, fractures_fr, b_fractures]

    print(fracture_mesh_step)
    #fractures_fr.set_mesh_step(fracture_mesh_step)

    factory.keep_only(*mesh_groups)
    factory.remove_duplicate_entities()
    factory.write_brep()

    min_el_size = fracture_mesh_step / 10
    fracture_el_size = np.max(dimensions) / 20
    max_el_size = np.max(dimensions) / 8


    fracture_el_size = field.constant(fracture_mesh_step, 10000)
    frac_el_size_only = field.restrict(fracture_el_size, fractures_fr, add_boundary=True)
    field.set_mesh_step_field(frac_el_size_only)

    mesh = options.Mesh()
    #mesh.Algorithm = options.Algorithm2d.MeshAdapt # produce some degenerated 2d elements on fracture boundaries ??
    #mesh.Algorithm = options.Algorithm2d.Delaunay
    #mesh.Algorithm = options.Algorithm2d.FrontalDelaunay
    #mesh.Algorithm3D = options.Algorithm3d.Frontal
    #mesh.Algorithm3D = options.Algorithm3d.Delaunay
    mesh.ToleranceInitialDelaunay = 0.01
    #mesh.ToleranceEdgeLength = fracture_mesh_step / 5
    mesh.CharacteristicLengthFromPoints = True
    mesh.CharacteristicLengthFromCurvature = True
    mesh.CharacteristicLengthExtendFromBoundary = 2
    mesh.CharacteristicLengthMin = min_el_size
    mesh.CharacteristicLengthMax = max_el_size
    mesh.MinimumCirclePoints = 6
    mesh.MinimumCurvePoints = 2


    #factory.make_mesh(mesh_groups, dim=2)
    factory.make_mesh(mesh_groups)
    factory.write_mesh(format=gmsh.MeshFormat.msh2)
    os.rename(mesh_name + ".msh2", mesh_file)
    #factory.show()