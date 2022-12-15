import logging

import numpy as np
import pandas as pd
from typing import *
import matplotlib
import matplotlib.pyplot as plt
import pyvista as pv
from matplotlib import ticker, cm

from .common import File
from endorse.indicator import IndicatorFn, Indicator
from mlmc.moments import Legendre
from mlmc.estimator import Estimate
from mlmc.quantity.quantity import make_root_quantity
from mlmc.sample_storage_hdf import SampleStorageHDF
from mlmc.quantity.quantity_estimate import estimate_mean


def plot_field(points, values, cut=(1,2), file=None):
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.set_aspect('equal', 'box')
    #ax.set_ylim(-50, 50)
    #ax.set_xlim(-50, 50)

    axis_label=('x','y','z')
    ax.set_xlabel(f"${axis_label[cut[0]]}$", fontsize=20)
    ax.set_ylabel(f"${axis_label[cut[1]]}$", fontsize=20)
    assert points.shape[0] == len(values)
    if len(values) > 1000:
        subset = np.random.randint(0,len(values),size=1000)
        points = points[subset, :]
        values = values[subset]
    sc = ax.scatter(points[:, cut[0]], points[:, cut[1]], c=values, s=1, cmap=plt.cm.viridis)

    axx = ax.twinx()
    isort = np.argsort(points[:,cut[0]])
    X = points[isort,cut[0]]
    Y = values[isort]
    axx.plot(X, Y)
    # levels = np.array([])
    #c = ax.contourf(X, Y, porosity, cmap=plt.cm.viridis)
    cb = fig.colorbar(sc)
    if file:
        fig.savefig(file)
    else:
        plt.show()


def plot_source(source_file):
    df = pd.read_csv(source_file.path)
    fig, ax = plt.subplots()
    times = np.array(df.iloc[:, 0])
    conc_flux = np.array(df.iloc[:, 1])
    ax.plot(times, conc_flux)
    ax.set_xlabel(df.columns[0])
    ax.set_ylabel(df.columns[1], color='blue')
    ax.set_xscale('log')
    ax.set_yscale('log')

    mass = (conc_flux[1:] + conc_flux[:-1])/2 * (times[1:] - times[:-1])
    cumul_mass = np.cumsum(mass)

    ax1 = ax.twinx()
    ax1.plot(times[1:], cumul_mass, c='red')
    ax1.set_ylabel('mass [kg]', color='red')
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-1, 1))
    ax1.yaxis.set_major_formatter(formatter)
    fig.tight_layout()
    fig.savefig("source_plot.pdf")
    #plt.show()


def plot_indicators(ind_functions: List[IndicatorFn], file=None):
    matplotlib.rcParams.update({'font.size': 16})
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']

    fig, ax = plt.subplots(figsize=(8, 6))
    for i, ind_fn in enumerate(ind_functions):
        tmax, vmax = ind_fn.time_max()

        #label = f"{ind_fn.indicator.indicator_label}; max: ({vmax:.2e}, {tmax:.2e})"
        label = f"{ind_fn.indicator.indicator_label}"
        ax.plot(ind_fn.times_fine(), ind_fn.spline(ind_fn.times_fine()), c=colors[i], label=label)
        ax.scatter(ind_fn.times, ind_fn.ind_values, marker='.', c=colors[i])
        ax.scatter([tmax], [vmax], s=100, c=colors[i], marker='*')
        #plt.text(tmax, vmax, f'({tmax:.2e}, {vmax:.2e})')
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-1, 1))
    ax.yaxis.set_major_formatter(formatter)

    yf = ticker.ScalarFormatter(useMathText=True)
    yf.set_scientific(True)
    yf.set_powerlimits((3, 3))
    ax.xaxis.set_major_formatter(yf)

    ax.set_xlabel("years")
    ax.set_ylabel("conc [g/m3]")
    plt.legend(loc='best')
    fig.tight_layout()
    if file is None:
        file = "indicators_plot.pdf"
    else:
        file = f"{file}.pdf"
    fig.savefig(file)


def _get_samples(quantity, sample_storage):
    n_moments = 5
    estimated_domain = Estimate.estimate_domain(quantity, sample_storage, quantile=0.001)
    moments_fn = Legendre(n_moments, estimated_domain)
    estimator = Estimate(quantity=quantity, sample_storage=sample_storage, moments_fn=moments_fn)
    samples = estimator.get_level_samples(level_id=0)[..., 0]
    return samples


def _get_values(hdf5_path):
    sample_storage = SampleStorageHDF(file_path=hdf5_path)
    sample_storage.chunk_size = 1024
    result_format = sample_storage.load_result_format()
    root_quantity = make_root_quantity(sample_storage, result_format)

    conductivity = root_quantity['indicator_conc']
    time = conductivity[1]  # times: [1]
    location = time['0']  # locations: ['0']
    values = location  # result shape: (10, 1)
    values = values[:4]
    values = values.select(values < 1e-1)
    samples = _get_samples(values, sample_storage)

    values = np.log(values)
    q_mean = estimate_mean(values)
    val_squares = estimate_mean(np.power(values - q_mean.mean, 2))
    std = np.sqrt(val_squares.mean)

    return q_mean.mean, std, samples


def plot_quantile_errorbar(data_dict, quantiles):
    matplotlib.rcParams.update({'font.size': 22})

    fig, ax = plt.subplots(figsize=(12, 10))

    pos = np.array([0, 0.1, 0.2, 0.3])
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']
    xticks_pos = []

    err_bars = []

    for case, hdf5_path in data_dict.items():
        mean, std, samples = _get_values(hdf5_path)
        exp_mean = np.exp(mean)
        yerr = np.exp(mean + np.array([-std, +std])) - exp_mean
        for i in [1,3]:
            # add samples
            s_pos = np.ones(len(samples[i])) * pos[i]
            ax.scatter(s_pos, samples[i], color=colors[i], marker="v")
            ax.set_yscale('log')

            # add errorbar, not able to pass array of colors for every quantile case
            err = ax.errorbar(pos[i], exp_mean[i], yerr=([-yerr[0, i]], [+yerr[1, i]]),
                      lw=2, capsize=6, capthick=2,
                      c=colors[i], marker="o", markersize=8, fmt=' ', linestyle='')

        xticks_pos.append(pos[int(len(pos)/2)])  # works sufficiently for x labels' centering

        pos += 1

    ax.set_ylabel('conc ' + r'$[g/m^3]$')
    ax.set_xticks(xticks_pos)
    ax.set_xticklabels(list(data_dict.keys()))

    labels = []
    for i in [1,3]:
        q = quantiles[i]
        exp = int(np.floor(np.log10(q)))
        man = int(q / 10 ** exp)
        label = f"$1 - {man}\\times 10^{{{exp:1d}}}$"

        #ind = Indicator.quantile(q_i)
        labels.append(label)

    ax.legend(labels=labels, loc=1)

    #ax.set_yscale("log")
    plt.savefig("quantiles.pdf")
    #plt.show()


def plot_contour(polydata: List[Any], attr_name: str, subtitles: List[str] = [None], title: str = None, file=None):
    nplots = len(polydata)
    # annotations of each plot
    use_titles = True
    if len(subtitles) != nplots:
        use_titles = False
        print('Skipping titles, mismatch with data length')
    # normalization
    vmin = min(min(data.point_data[f'{attr_name}']) for data in polydata)
    vmax = max(max(data.point_data[f'{attr_name}']) for data in polydata)
    # setup figure
    fig, ax = plt.subplots(nplots, 1, figsize=(8, 6))
    if title:
        fig.suptitle(title, fontsize=14)
    cs = [None] * nplots
    for i, data in enumerate(polydata):
        x = data.points
        tri = data.faces.reshape((-1, 4))[:, 1:]
        u = data.point_data[f'{attr_name}']
        # shift axis to readable coordinates and display origin
        xaxshift, yaxshift = min(x[:, 0]), min(x[:, 1])
        # plot contour with common colormap
        cs[i] = ax[i].tricontourf(x[:, 0] - xaxshift, x[:, 1] - yaxshift, tri, u, vmin=vmin, vmax=vmax,
                                  cmap=cm.coolwarm)
        #     unorm = (u-vmin+1e-10)/max(u-vmax)
        #     ax.tricontourf(x[:,0], x[:,1], tri, unorm, locator=ticker.LogLocator())
        # setup annotations
        formatter = ticker.ScalarFormatter(useMathText=True)
        ax[i].yaxis.set_major_formatter(formatter)
        ax[i].set_ylabel("y[-]", fontsize=10)
        ax[i].set_xticks([])

        if use_titles:
            ax[i].legend([subtitles[i]], fontsize=6)

    ax[-1].xaxis.set_major_locator(ticker.AutoLocator())
    ax[-1].set_xlabel("x[-]", fontsize=10)
    ax[-1].text(-3, -5, f'[{xaxshift:.1e}, \n{yaxshift:.1e}]', ha='right', va='top', fontsize=10)
    # add colorbar
    fig.subplots_adjust(right=0.825, hspace=0.1, top=0.92)
    cax = fig.add_axes([ax[-1].get_position().x1 + 0.02, ax[-1].get_position().y0, 0.035,
                        ax[0].get_position().y1 - ax[-1].get_position().y0])
    m = cm.ScalarMappable(cmap=cm.coolwarm)
    m.set_clim(vmin, vmax)
    cbar = fig.colorbar(m, cax=cax)
    cbar.set_label('conc [g/m3]', rotation=270, fontsize=10, labelpad=5)
    # file output
    if file is None:
        file = "slices_plot.pdf"
    else:
        file = f"{file}.pdf"
    fig.savefig(file)


def plot_slices(pvd_in: File, attr_name: str, z_loc, timeidxs):
    pvd_content = pv.get_reader(pvd_in.path)
    times = np.asarray(pvd_content.time_values)
    for i in timeidxs:
        if not 0 <= i < len(times):
            logging.warning(f'Cannot plot slice, index {i} out of range')
            continue
        pvd_content.set_active_time_point(i)
        dataset = pvd_content.read()
        title = f'Concentration at t={times[i]:.1e}'
        plane = []
        subtitles = []
        for z in z_loc:
            # assumption of one block for polyblock
            p = dataset.slice(normal=[0, 0, 1], origin=[0, 0, z], generate_triangles=True)[0]
            plane.append(p)
            subt = f'z={z}'
            subtitles.append(subt)
        plot_contour(plane, attr_name, subtitles=subtitles, title=title)
