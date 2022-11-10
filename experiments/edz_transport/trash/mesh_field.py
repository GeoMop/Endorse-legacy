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

#def apply_field(field1, field2, reference_fn, dim=2, tolerance=0.15, max_mismatch=5, mesh_name="field_mesh"):
def apply_field(dim, tolerance=0.15, max_mismatch=5, mesh_name="field_mesh"):
    """
    Create a mesh of dimension dim on a unit cube and
    compare element sizes to given reference function of coordinates.
    """

    sample_dir = "output"
    os.chdir(sample_dir)

    model = gmsh.GeometryOCC(mesh_name)

    #circ = model.circle(10).set_region("circle")

    rec3 = model.rectangle([220, 20]).set_region("square3")
    rec = model.rectangle([300, 100]).cut(rec3.copy()).set_region("square")
    #rec = model.rectangle([500, 500]).cut(rec3.copy(), rec4.copy()).set_region("square")
    #rec = model.group(rec3, rec2)

    boundaries = gmsh.ObjectSet.get_boundary(rec3, False)
    boundaries2 = gmsh.ObjectSet.get_boundary(rec, False)

    #rec3 = model.rectangle([50, 10]).translate([100,50,0])

    #tagy = rec1.tags

    a = 3
    b = 0.2
    c = 10
    d = 10

    field3 = field.Field('MathEval')
    # field3.F = "x/100"

    uzly = []
    for i in range(boundaries.size):
        uzly.append(boundaries.dim_tags[i][1])

    #for j in range(boundaries2.size):
    #   uzly.append(boundaries2.dim_tags[i][1])


    #field1 = field.threshold(field.distance_nodes([boundaries.dim_tags[0][1], boundaries.dim_tags[1][1]]), lower_bound=(a, b), upper_bound=(c, d))
    f = field.threshold(field.distance_edges(uzly, 30), lower_bound=(a, b), upper_bound=(c, d))

    """
    pta = [0, 0, 0]
    ptb = [100, 100, 0]

    field2 = field.box(pta, ptb, 1, 20)

    f_1 = field.constant(1)
    #f1 = (field.abs(field.y)- 10)/10
    f1 = (field.abs(field.y) - 10)
    f2 = (field.abs(field.x) - 50)/10

    f1 = (1 - field.sign(field.abs(field.y) - 10)) + (1 + field.sign(field.abs(field.y) - 10)) * field.abs(field.abs(field.y) - 10)/8
    f2 = (1 - field.sign(field.abs(field.x) - 50)) + (1 + field.sign(field.abs(field.x) - 50)) * field.abs(
        field.abs(field.x) - 50) / 8

    f1 = 5 - field.y * field.y / 100 - field.x * field.x / 100
    f2 = field.y * field.y / 500 + field.x * field.x / 500

    #f2 = (field.abs(field.y) +10) / 10

    f = field.maximum(f1, f2)
    f = field1
    """

    ####################################################
    """
    model = gmsh.GeometryOCC(mesh_name)

    rec3 = model.box([200, 200, 200]).set_region("square3")
    rec4 = model.box([50, 50, 50]).set_region("square4")
    rec = rec3.cut(rec4.copy()).set_region("square")

    boundaries = gmsh.ObjectSet.get_boundary(rec4, False)

    uzly = []
    for i in range(boundaries.size):
        uzly.append(boundaries.dim_tags[i][1])

    f = field.threshold(field.distance_surfaces(uzly, 8), lower_bound=(a, b), upper_bound=(c, d))
    """
    ####################################################

    model.set_mesh_step_field(f)

    model.write_brep()
    model.mesh_options.CharacteristicLengthMin = 0.1
    model.mesh_options.CharacteristicLengthMax = 100
    model.make_mesh([rec], dim=dim)
    model.write_mesh(mesh_name + ".msh2", gmsh.MeshFormat.msh2)

    del model
