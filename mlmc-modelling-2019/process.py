import sys
import os
script_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(script_dir, '../MLMC/src'))
sys.path.append(os.path.join(script_dir, '../gmsh-tools/src'))
sys.path.append(os.path.join(script_dir, '../PyBS/src'))

import shutil
import subprocess
import yaml
import attr
import numpy as np
import collections
import gmsh_io
# import matplotlib.pyplot as plt

import mlmc.random.fracture as fracture

# TODO:
# - enforce creation of empty physical groups, or creation of empty regions in the flow input
# - speedup mechanics

@attr.s(auto_attribs=True)
class ValueDescription:
    time: float
    position: str
    quantity: str
    unit: str



def substitute_placeholders(file_in, file_out, params):
    """
    Substitute for placeholders of format '<name>' from the dict 'params'.
    :param file_in: Template file.
    :param file_out: Values substituted.
    :param params: { 'name': value, ...}
    """
    used_params = []
    with open(file_in, 'r') as src:
        text = src.read()
    for name, value in params.items():
        placeholder = '<%s>' % name
        n_repl = text.count(placeholder)
        if n_repl > 0:
            used_params.append(name)
            text = text.replace(placeholder, str(value))
    with open(file_out, 'w') as dst:
        dst.write(text)
    return used_params

def to_polar(x, y, z):
    rho = np.sqrt(x ** 2 + y ** 2)
    phi = np.arctan2(y, x)
    if z > 0:
        phi += np.pi
    return (phi, rho)

def plot_fr_orientation(fractures):
    family_dict = collections.defaultdict(list)
    for fr in fractures:
        x, y, z = fracture.FisherOrientation.rotate(np.array([0,0,1]), axis=fr.rotation_axis, angle=fr.rotation_angle)[0]
        family_dict[fr.region].append([
            to_polar(z, y, x),
            to_polar(z, x, -y),
            to_polar(y, x, z)
            ])

    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, subplot_kw=dict(projection='polar'))
    for name, data in family_dict.items():
        # data shape = (N, 3, 2)
        data = np.array(data)
        for i, ax in enumerate(axes):
            phi = data[:, i, 0]
            r = data[:, i, 1]
            c = ax.scatter(phi, r, cmap='hsv', alpha=0.75, label=name)
    axes[0].set_title("X-view, Z-north")
    axes[1].set_title("Y-view, Z-north")
    axes[2].set_title("Z-view, Y-north")
    for ax in axes:
        ax.set_theta_zero_location("N")
        ax.set_theta_direction(-1)
        ax.set_ylim(0, 1)
    fig.legend(loc = 1)
    fig.savefig("fracture_orientation.pdf")
    plt.close(fig)
    #plt.show()

def generate_fractures(config_dict):
    geom = config_dict["geometry"]
    lx, ly = geom["fractures_box"]
    fr_size_range = geom["pow_law_size_range"]
    pow_law_exp_3d = geom["pow_law_size_exp"]
    pow_law_sample_range = geom["pow_law_sample_range"]
    n_frac_limit = geom["n_frac_limit"]
    p_32 = geom["p_32"]

    # generate fracture set
    fracture_box = [lx, ly, 0]
    area = lx * ly

    pop = fracture.Population(area, fracture.LineShape)
    pop.add_family("all",
                   fracture.FisherOrientation(0, 90, np.inf),
                   fracture.VonMisesOrientation(0, 0),
                   fracture.PowerLawSize.from_mean_area(pow_law_exp_3d-1, fr_size_range, p_32, pow_law_exp_3d))

    if pow_law_sample_range:
        pop.set_sample_range(pow_law_sample_range)
    elif n_frac_limit:
        pop.set_sample_range([None, max(lx, ly)], sample_size=n_frac_limit)

    print("total mean size: ", pop.mean_size())
    print("size range:", pop.families[0].size.sample_range)
    pos_gen = fracture.UniformBoxPosition(fracture_box)
    fractures = pop.sample(pos_distr=pos_gen, keep_nonempty=True)


    # for fr in fractures:
    #     fr.region = "fr"
    # used_families = set((f.region for f in fractures))
    # for model in ["hm_params", "th_params", "th_params_ref"]:
    #     model_dict = config_dict[model]
    #     model_dict["fracture_regions"] = list(used_families)
    #     model_dict["left_well_fracture_regions"] = [".{}_left_well".format(f) for f in used_families]
    #     model_dict["right_well_fracture_regions"] = [".{}_right_well".format(f) for f in used_families]
    fr_set = fracture.Fractures(fractures, fr_size_range[0] / 2)
    #fr_set.fragment()
    #fr_set.simplify()
    return fr_set




def fractures_from_base_shape(gmsh_geom, fractures, base_shape: 'ObjectSet'):
    # Copy and transform the given `base_shape` into set of fracture shapes
    # according to the fracture properties given by `fractures` list.
    # Return list of fragmented shapes.
    shapes = []
    for i, fr in enumerate(fractures):
        shape = base_shape.copy()
        print("fr: ", i, "tag: ", shape.dim_tags, "axis:", fr.rotation_axis, "angle:", fr.rotation_angle)
        shape = shape.scale([fr.rx, fr.ry, 1]) \
            .rotate(axis=[0, 0, 1], angle=fr.shape_angle) \
            .rotate(axis=fr.rotation_axis, angle=fr.rotation_angle) \
            .translate(fr.centre) \
            .set_region(fr.region)

        shapes.append(shape)

    fracture_fragments = gmsh_geom.fragment(*shapes)
    return fracture_fragments

def fractures_lines(gmsh_geom, fractures):
    # Copy and transform the given `base_shape` into set of fracture shapes
    # according to the fracture properties given by `fractures` list.
    # Return list of fragmented shapes.
    shapes = []
    base_line = np.array([[-0.5, 0, 0], [0.5, 0, 0]])
    for i, fr in enumerate(fractures):
        line = fracture.FisherOrientation.rotate(base_line * fr.rx, np.array([0,0,1]), fr.shape_angle)
        line += fr.centre
        shape = gmsh_geom.line(line[0], line[1]).set_region(fr.region)
        shapes.append(shape)
        print("fr: ", i, "tag: ", shape.dim_tags)


    fracture_fragments = gmsh_geom.fragment(*shapes)
    return fracture_fragments


def create_fractures_polygons(gmsh_geom, fractures):
    # From given fracture date list 'fractures'.
    # transform the base_shape to fracture objects
    # fragment fractures by their intersections
    # return dict: fracture.region -> GMSHobject with corresponding fracture fragments
    frac_obj = fracture.Fractures(fractures)
    frac_obj.snap_vertices_and_edges()
    shapes = []
    for fr, square in zip(fractures, frac_obj.squares):
        shape = gmsh_geom.make_polygon(square).set_region(fr.region)
        shapes.append(shape)

    fracture_fragments = gmsh_geom.fragment(*shapes)
    return fracture_fragments


def make_dfn_mesh(config_dict, fractures, mesh_name, mesh_file):
    # 2d Mesh containing just the fractures, no bulk.

    geom = config_dict["geometry"]
    fracture_mesh_step = geom['fracture_mesh_step']
    lx, ly = geom["domain_box"]

    print("load gmsh api")

    from gmsh_api import gmsh
    from gmsh_api import options
    from gmsh_api import field

    factory = gmsh.GeometryOCC(mesh_name, verbose=True)
    gopt = options.Geometry()
    gopt.Tolerance = 0.0001
    gopt.ToleranceBoolean = 0.001
    # gopt.MatchMeshTolerance = 1e-1

    # Main box

    domain = factory.rectangle([lx, ly])
    sides = dict(
        side_y0=factory.line([-lx/2, -ly/2, 0], [+lx/2, -ly/2, 0]),
        side_y1=factory.line([-lx/2, +ly/2, 0], [+lx/2, +ly/2, 0]),
        side_x0=factory.line([-lx/2, -ly/2, 0], [-lx/2, +ly/2, 0]),
        side_x1=factory.line([+lx/2, -ly/2, 0], [+lx/2, +ly/2, 0])
    )
    for name, side in sides.items():
        side.modify_regions(name)

    print("n fractures:", len(fractures))
    base_line_shape = factory.line([-0.5, 0, 0], [0.5, 0, 0])
    #fractures = fractures_from_base_shape(factory, fractures, base_line_shape)
    fractures = fractures_lines(factory, fractures)
    fractures_group = factory.group(*fractures)
    #fractures_group = fractures_group.remove_small_mass(fracture_mesh_step * fracture_mesh_step / 10)


    # fractures, fragmented, fractures boundary
    # print("cut fractures by the domain")
    # fractures_group = fractures_group.intersect(domain.copy())
    # print("select boundaries")
    # b_domain_fr = fractures_group.get_boundary()
    mesh_objects = [fractures_group]
    # for name, side_tool in sides.items():
    #     isec = b_domain_fr.select_by_intersect(side_tool)
    #     mesh_objects.append(isec.modify_regions("." + name))


    print(fracture_mesh_step)
    #fractures_fr.set_mesh_step(fracture_mesh_step)

    factory.keep_only(*mesh_objects)
    print("remove duplicities")
    factory.remove_duplicate_entities()
    factory.write_brep()

    min_el_size = fracture_mesh_step / 10
    max_el_size = np.max([lx, ly]) / 8

    fracture_el_size = field.constant(fracture_mesh_step, 10000)
    frac_el_size_only = field.restrict(fracture_el_size, fractures_group, add_boundary=False)
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
    factory.make_mesh(mesh_objects, dim=1)
    factory.write_mesh(format=gmsh.MeshFormat.msh2)
    os.rename(mesh_name + ".msh2", mesh_file)
    factory.show()



def make_2d_1d_mesh(config_dict, fractures, mesh_name, mesh_file):
    # 2d Mesh containing just the fractures, no bulk.

    geom = config_dict["geometry"]
    fracture_mesh_step = geom['fracture_mesh_step']
    lx, ly = geom["domain_box"]

    print("load gmsh api")

    from gmsh_api import gmsh
    from gmsh_api import options
    from gmsh_api import field

    factory = gmsh.GeometryOCC(mesh_name, verbose=True)
    gopt = options.Geometry()
    gopt.Tolerance = 0.0001
    gopt.ToleranceBoolean = 0.001
    # gopt.MatchMeshTolerance = 1e-1

    # Main box

    domain = factory.rectangle([lx, ly])
    sides = dict(
        side_y0=factory.line([-lx/2, -ly/2, 0], [+lx/2, -ly/2, 0]),
        side_y1=factory.line([-lx/2, +ly/2, 0], [+lx/2, +ly/2, 0]),
        side_x0=factory.line([-lx/2, -ly/2, 0], [-lx/2, +ly/2, 0]),
        side_x1=factory.line([+lx/2, -ly/2, 0], [+lx/2, +ly/2, 0])
    )
    for name, side in sides.items():
        side.modify_regions(name)

    print("n fractures:", len(fractures))
    base_line_shape = factory.line([-0.5, 0, 0], [0.5, 0, 0])
    #fractures = fractures_from_base_shape(factory, fractures, base_line_shape)
    fractures = fractures_lines(factory, fractures)
    fractures_group = factory.group(*fractures)
    #fractures_group = fractures_group.remove_small_mass(fracture_mesh_step * fracture_mesh_step / 10)


    # fractures, fragmented, fractures boundary
    print("cut fractures by the domain")
    fractures_group = fractures_group.intersect(domain.copy())
    domain_fr, fractures_fr = factory.fragment(domain, fractures_group)

    print("select boundaries")
    b_fr = fractures_fr.get_boundary()
    b_domain = domain_fr.get_boundary()
    mesh_objects = [domain_fr, fractures_fr]
    for side_name, side_tool in sides.items():
        isec = b_fr.select_by_intersect(side_tool)
        mesh_objects.append(isec.modify_regions(".{}_" + side_name))
        isec = b_domain.select_by_intersect(side_tool)
        mesh_objects.append(isec.modify_regions(".{}_" + side_name))

    print(fracture_mesh_step)
    fractures_fr.set_mesh_step(fracture_mesh_step)
    #b_domain.set_mesh_step(fracture_mesh_step)

    factory.keep_only(*mesh_objects)
    print("remove duplicities")
    factory.remove_duplicate_entities()
    factory.write_brep()

    min_el_size = fracture_mesh_step
    max_el_size = np.max([lx, ly]) / 8

    #fracture_el_size = field.constant(fracture_mesh_step, 10000)
    #frac_el_size_only = field.restrict(fracture_el_size, fractures_fr, add_boundary=False)
    #field.set_mesh_step_field(frac_el_size_only)

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


    factory.make_mesh(mesh_objects, dim=2)
    factory.write_mesh(format=gmsh.MeshFormat.msh2)
    os.rename(mesh_name + ".msh2", mesh_file)
    factory.show()



def make_dfn_mesh_direct(geom, fractures, mesh_base):
    import frac_geom as fg
    geom = config_dict["geometry"]

    fracture_mesh_step = geom['fracture_mesh_step']
    bulk_mesh_step = fracture_mesh_step
    lx, ly = geom["domain_box"]
    root_polygon = [[-lx / 2, -ly / 2], [+lx / 2, -ly / 2], [+lx / 2, +ly / 2], [-lx / 2, +ly / 2]]
    frac_lines = np.array([(fractures.points[p0][:2], fractures.points[p1][:2]) for p0, p1 in fractures.lines])
    mesh = fg.make_frac_mesh(root_polygon, bulk_mesh_step, frac_lines, fracture_mesh_step, mesh_base=mesh_base)
    return mesh

"""
TODO:
- check if the mesh file is created
- return the mesh file
- try with moderate num of fracs
- try with higher number of frac
- compute fracture conductivities
- run flow
"""

def prepare_mesh(config_dict, fractures):
    mesh_name = config_dict["mesh_name"]
    mesh_file = mesh_name + ".msh"
    if True: #not os.path.isfile(mesh_file):
        #make_dfn_mesh(config_dict, fractures, mesh_name, mesh_file)
        #make_2d_1d_mesh(config_dict, fractures, mesh_name, mesh_file)
        healed_mesh = make_dfn_mesh_direct(config_dict, fractures, mesh_name)
    mesh_healed = mesh_name + "_healed.msh"
    # if not os.path.isfile(mesh_healed):
    #     import heal_mesh
    #     hm = heal_mesh.HealMesh.read_mesh(mesh_file, node_tol=1e-4)
    #     hm.heal_mesh(gamma_tol=0.01)
    #     hm.stats_to_yaml(mesh_name + "_heal_stats.yaml")
    #     hm.write()
    #     assert hm.healed_mesh_name == mesh_healed
    return mesh_healed

def check_conv_reasons(log_fname):
    with open(log_fname, "r") as f:
        for line in f:
            tokens = line.split(" ")
            try:
                i = tokens.index('convergence')
                if tokens[i + 1] == 'reason':
                    value = tokens[i + 2].rstrip(",")
                    conv_reason = int(value)
                    if conv_reason < 0:
                        print("Failed to converge: ", conv_reason)
                        return False
            except ValueError:
                continue
    return True

def call_flow(config_dict, param_key, result_files):
    """
    Redirect sstdout and sterr, return true on succesfull run.
    :param arguments:
    :return:
    """

    params = config_dict[param_key]
    fname = params["in_file"]
    substitute_placeholders(fname + '_tmpl.yaml', fname + '.yaml', params)
    arguments = config_dict["_aux_flow_path"].copy()
    output_dir = "output_" + fname
    config_dict[param_key]["output_dir"] = output_dir
    if all([os.path.isfile(os.path.join(output_dir, f)) for f in result_files]):
        status = True
    else:
        arguments.extend(['--output_dir', output_dir, fname + ".yaml"])
        print("Running: ", " ".join(arguments))
        with open(fname + "_stdout", "w") as stdout:
            with open(fname + "_stderr", "w") as stderr:
                completed = subprocess.run(arguments, stdout=stdout, stderr=stderr)
        print("Exit status: ", completed.returncode)
        status = completed.returncode == 0
    conv_check = check_conv_reasons(os.path.join(output_dir, "flow123.0.log"))
    print("converged: ", conv_check)
    return status # and conv_check



def prepare_th_input(config_dict):
    """
    Prepare FieldFE input file for the TH simulation.
    :param config_dict: Parsed config.yaml. see key comments there.
    """
    # pass
    # we have to read region names from the input mesh
    # input_mesh = gmsh_io.GmshIO(config_dict['hm_params']['mesh'])
    #
    # is_bc_region = {}
    # for name, (id, _) in input_mesh.physical.items():
    #     unquoted_name = name.strip("\"'")
    #     is_bc_region[id] = (unquoted_name[0] == '.')

    # read mesh and mechanichal output data
    mechanics_output = os.path.join(config_dict['hm_params']["output_dir"], 'mechanics.msh')
    mesh = gmsh_io.GmshIO(mechanics_output)

    n_bulk = len(mesh.elements)
    ele_ids = np.array(list(mesh.elements.keys()), dtype=float)

    init_fr_cs = float(config_dict['hm_params']['fr_cross_section'])
    init_fr_K = float(config_dict['hm_params']['fr_conductivity'])
    init_bulk_K = float(config_dict['hm_params']['bulk_conductivity'])
    min_fr_cross_section = float(config_dict['th_params']['min_fr_cross_section'])
    max_fr_cross_section = float(config_dict['th_params']['max_fr_cross_section'])

    time_idx = 1
    time, field_cs = mesh.element_data['cross_section_updated'][time_idx]

    # cut small and large values of cross-section
    cs = np.maximum(np.array([v[0] for v in field_cs.values()]), min_fr_cross_section)
    cs = np.minimum(cs, max_fr_cross_section)

    K = np.where(
        cs == 1.0,      # condition
        init_bulk_K,    # true array
        init_fr_K * (cs / init_fr_cs) ** 2
    )

    # get cs and K on fracture elements only
    fr_indices = np.array([int(key) for key, val in field_cs.items() if val[0] != 1])
    cs_fr = np.array([cs[i] for i in fr_indices])
    k_fr = np.array([K[i] for i in fr_indices])

    # compute cs and K statistics and write it to a file
    fr_param = {}
    avg = float(np.average(cs_fr))
    median = float(np.median(cs_fr))
    interquantile = float(1.5 * (np.quantile(cs_fr, 0.75) - np.quantile(cs_fr, 0.25)))
    fr_param["fr_cross_section"] = {"avg": avg, "median": median, "interquantile": interquantile}

    avg = float(np.average(k_fr))
    median = float(np.median(k_fr))
    interquantile = float(1.5 * (np.quantile(k_fr, 0.75) - np.quantile(k_fr, 0.25)))
    fr_param["fr_conductivity"] = {"avg": avg, "median": median, "interquantile": interquantile}

    with open('fr_param_output.yaml', 'w') as outfile:
        yaml.dump(fr_param, outfile, default_flow_style=False)

    # mesh.write_fields('output_hm/th_input.msh', ele_ids, {'conductivity': K})
    th_input_file = 'th_input.msh'
    with open(th_input_file, "w") as fout:
        mesh.write_ascii(fout)
        mesh.write_element_data(fout, ele_ids, 'conductivity', K[:, None])
        mesh.write_element_data(fout, ele_ids, 'cross_section_updated', cs[:, None])

    # create field for K (copy cs)
    # posun dat K do casu 0
    # read original K = oK (define in config yaml)
    # read original cs = ocs (define in config yaml)
    # compute K = oK * (cs/ocs)^2
    # write K

    # posun dat cs do casu 0
    # write cs

    # mesh.element_data.



# def get_result_description():
#     """
#     :return:
#     """
#     end_time = 30
#     values = [ [ValueDescription(time=t, position="extraction_well", quantity="power", unit="MW"),
#                 ValueDescription(time=t, position="extraction_well", quantity="temperature", unit="Celsius deg.")
#                 ] for t in np.linspace(0, end_time, 0.1)]
#     power_series, temp_series = zip(*values)
#     return power_series + temp_series
#
#
# def extract_time_series(yaml_stream, regions, extract):
#     """
#
#     :param yaml_stream:
#     :param regions:
#     :return: times list, list: for every region the array of value series
#     """
#     data = yaml.safe_load(yaml_stream)['data']
#     times = set()
#     reg_series = {reg: [] for reg in regions}
#
#     for time_data in data:
#         region = time_data['region']
#         if region in reg_series:
#             times.add(time_data['time'])
#             power_in_time = extract(time_data)
#             reg_series[region].append(power_in_time)
#     times = list(times)
#     times.sort()
#     series = [np.array(region_series, dtype=float) for region_series in reg_series.values()]
#     return np.array(times), series
#
#
# def extract_results(config_dict):
#     """
#     :param config_dict: Parsed config.yaml. see key comments there.
#     : return
#     """
#     bc_regions = ['.fr_left_well', '.left_well', '.fr_right_well', '.right_well']
#     out_regions = bc_regions[2:]
#     output_dir = config_dict["th_params"]["output_dir"]
#     with open(os.path.join(output_dir, "energy_balance.yaml"), "r") as f:
#         power_times, reg_powers = extract_time_series(f, bc_regions, extract=lambda frame: frame['data'][0])
#         power_series = -sum(reg_powers)
#
#     with open(os.path.join(output_dir, "Heat_AdvectionDiffusion_region_stat.yaml"), "r") as f:
#         temp_times, reg_temps = extract_time_series(f, out_regions, extract=lambda frame: frame['average'][0])
#     with open(os.path.join(output_dir, "water_balance.yaml"), "r") as f:
#         flux_times, reg_fluxes = extract_time_series(f, out_regions, extract=lambda frame: frame['data'][0])
#     sum_flux = sum(reg_fluxes)
#     avg_temp = sum([temp * flux for temp, flux in zip(reg_temps, reg_fluxes)]) / sum_flux
#     print("temp: ", avg_temp)
#     return temp_times, avg_temp, power_times, power_series
#
#
# def plot_exchanger_evolution(temp_times, avg_temp, power_times, power_series):
#     abs_zero_temp = 273.15
#     year_sec = 60 * 60 * 24 * 365
#
#     import matplotlib.pyplot as plt
#     fig, ax1 = plt.subplots()
#     temp_color = 'red'
#     ax1.set_xlabel('time [y]')
#     ax1.set_ylabel('Temperature [C deg]', color=temp_color)
#     ax1.plot(temp_times[1:] / year_sec, avg_temp[1:] - abs_zero_temp, color=temp_color)
#     ax1.tick_params(axis='y', labelcolor=temp_color)
#
#     ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
#     pow_color = 'blue'
#     ax2.set_ylabel('Power [MW]', color=pow_color)  # we already handled the x-label with ax1
#     ax2.plot(power_times[1:] / year_sec, power_series[1:] / 1e6, color=pow_color)
#     ax2.tick_params(axis='y', labelcolor=pow_color)
#
#     fig.tight_layout()  # otherwise the right y-label is slightly clipped
#     plt.show()

def sample_mesh_repository(mesh_repository):
    """
    
    """
    mesh_file = np.random.choice(os.listdir(mesh_repository))
    healed_mesh = "random_fractures_healed.msh"
    shutil.copyfile(os.path.join(mesh_repository, mesh_file), healed_mesh)
    heal_ref_report = { 'flow_stats': { 'bad_el_tol': 0.01, 'bad_elements': [], 'bins': [], 'hist': []},
                        'gamma_stats': {'bad_el_tol': 0.01, 'bad_elements': [], 'bins': [], 'hist': []}}
    with open("random_fractures_heal_stats.yaml", "w") as f:
        yaml.dump(heal_ref_report, f)
    return healed_mesh


def setup_dir(config_dict, clean=False):
    for f in config_dict["copy_files"]:
        shutil.copyfile(os.path.join(script_dir, f), os.path.join(".", f))
    flow_exec = config_dict["flow_executable"].copy()
    # if not os.path.isabs(flow_exec[0]):
    #     flow_exec[0] = os.path.join(script_dir, flow_exec[0])
    config_dict["_aux_flow_path"] = flow_exec

def sample(config_dict):

    setup_dir(config_dict, clean=True)
    #mesh_repo = config_dict.get('mesh_repository', None)
    fractures = generate_fractures(config_dict)
    healed_mesh = prepare_mesh(config_dict, fractures)

    # healed_mesh_bn = os.path.basename(healed_mesh)
    # config_dict["dfn_flow_params"]["mesh"] = healed_mesh_bn
    #
    # succeed = call_flow(config_dict, 'hm_params', result_files=["mechanics.msh"])
    print("Finished")


if __name__ == "__main__":
    sample_dir = sys.argv[1]
    with open(os.path.join(script_dir, "config.yaml"), "r") as f:
        config_dict = yaml.safe_load(f)

    os.chdir(sample_dir)
    np.random.seed(1)
    sample(config_dict)
