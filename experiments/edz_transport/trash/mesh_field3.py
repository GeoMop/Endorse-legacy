import pytest
from gmsh import model as gmsh_model
from bgem.gmsh import gmsh, field, options
import numpy as np
import os
import math
import yaml


def create_fractures_rectangles(gmsh_geom, fractures, base_shape: 'ObjectSet'):
    # From given fracture date list 'fractures'.
    # transform the base_shape to fracture objects
    # fragment fractures by their intersections
    # return dict: fracture.region -> GMSHobject with corresponding fracture fragments
    shapes = []
    for i, fr in enumerate(fractures):
        shape = base_shape.copy()
        print("fr: ", i, "tag: ", shape.dim_tags)
        shape = shape.scale([fr.rx, fr.ry, 1]) \
            .rotate(axis=fr.rotation_axis, angle=fr.rotation_angle) \
            .translate(fr.centre) \
            .set_region(fr.region)

        shapes.append(shape)

    fracture_fragments = gmsh_geom.fragment(*shapes)
    return fracture_fragments

def apply_field3(dim, tolerance=0.15, max_mismatch=5, mesh_name="field_mesh"):
    """
    Create a mesh of dimension dim on a unit cube and
    compare element sizes to given reference function of coordinates.
    """

    sample_dir = "output"
    os.chdir(sample_dir)

    model = gmsh.GeometryOCC(mesh_name)

    rec = model.rectangle([20, 20]).set_region("square2")
    boundaries = gmsh.ObjectSet.get_boundary(rec, True)

    f1 = 3 - field.y * field.y / 5 - field.x * field.x / 5
    f2 = field.y * field.y / 50 + field.x * field.x / 50


    #f2 = (field.abs(field.y) +10) / 10

    f = field.maximum(f1, f2)
    #f = f2

    model.set_mesh_step_field(f)
    # model.set_mesh_step_field(field2)
    # rec1.mesh_step(1)
    # rec.mesh_step(0.5)

    model.write_brep()
    model.mesh_options.CharacteristicLengthMin = 0.5
    model.mesh_options.CharacteristicLengthMax = 10
    model.make_mesh([rec], dim=dim)
    model.write_mesh(mesh_name + ".msh2", gmsh.MeshFormat.msh2)

    del model

