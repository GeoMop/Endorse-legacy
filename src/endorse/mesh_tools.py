import os
from bgem.gmsh import field, options, gmsh

def edz_refinement_field(cfg:"dotdict", factory:"GeometryOCC") -> field.Field:
    """
    Refinement mesh step field for resolution of the EDZ.
    """
    b_cfg = cfg.borehole
    center_line = factory.line([0,0,0], [b_cfg.length, 0, 0]).translate([0, 0, b_cfg.z_pos])

    bx, by, bz = cfg.box_dimensions
    n_sampling = int(b_cfg.length / 2)
    dist = field.distance(center_line, sampling = n_sampling)
    inner = field.geometric(dist, a=(b_cfg.radius, cfg.edz_mesh_step * 0.9), b=(cfg.edz_radius, cfg.edz_mesh_step))
    outer = field.polynomial(dist, a=(cfg.edz_radius, cfg.edz_mesh_step), b=(by/2, cfg.boundary_mesh_step), q=1.7)
    return field.maximum(inner, outer)


def edz_meshing(cfg, factory, objects, mesh_file):
    """
    Common EDZ and transport domain meshing setup.
    """
    factory.write_brep()
    factory.mesh_options.CharacteristicLengthMin = cfg.get("min_mesh_step", cfg.boreholes_mesh_step)
    factory.mesh_options.CharacteristicLengthMax = cfg.boundary_mesh_step
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


def container_period(cfg):
    cont = cfg.containers
    return cont.length + cont.spacing


def container_x_pos(cfg, i_pos):
    cont = cfg.containers
    return cont.offset + i_pos * container_period(cfg)