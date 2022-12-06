import numpy as np
import pandas as pd
from typing import *
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import ticker

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
    plt.show()


def plot_indicators(ind_functions: List[IndicatorFn], file=None):
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']

    fig, ax = plt.subplots(figsize=(8, 6))
    for i, ind_fn in enumerate(ind_functions):
        tmax, vmax = ind_fn.time_max()

        label = f"{ind_fn.indicator.indicator_label}; max: ({vmax:.2e}, {tmax:.2e})"
        ax.plot(ind_fn.times_fine(), ind_fn.spline(ind_fn.times_fine()), c=colors[i], label=label)
        ax.scatter(ind_fn.times, ind_fn.ind_values, marker='.', c=colors[i])
        ax.scatter([tmax], [vmax], s=100, c=colors[i], marker='*')
        #plt.text(tmax, vmax, f'({tmax:.2e}, {vmax:.2e})')
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-1, 1))
    ax.yaxis.set_major_formatter(formatter)
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

    q_mean = estimate_mean(values)

    std = []
    for i in range(len(q_mean.mean)):
        std.append(np.sqrt(estimate_mean(np.power(values[i] - q_mean.mean[i], 2)).mean)[0])

    samples = _get_samples(values, sample_storage)

    return np.array(q_mean.mean[:4]).flatten(), std[:4], samples[:4]


def plot_quantile_errorbar(data_dict):
    matplotlib.rcParams.update({'font.size': 22})

    fig, ax = plt.subplots(figsize=(12, 10))

    pos = np.array([0, 0.1, 0.2, 0.3])
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']
    xticks_pos = []

    quantiles = [0.005, 0.002, 0.001, 0.0005]
    err_bars = []

    for case, hdf5_path in data_dict.items():
        mean, std, samples = _get_values(hdf5_path)
        for i in range(len(mean)):
            # add samples
            s_pos = np.ones(len(samples[i])) * pos[i]
            ax.scatter(s_pos, samples[i], color=colors[i], marker="v")

            # add errorbar
            err = ax.errorbar(pos[i], mean[i], yerr=std[i],
                              lw=2, capsize=6, capthick=2,
                              color=colors[i], marker="o", markersize=8)

            if len(err_bars) < len(pos):
                err_bars.append(err)

        xticks_pos.append(pos[int(len(pos)/2)])  # works sufficiently for x labels' centering
        pos += 1

    ax.set_ylabel('conc ' + r'$[g/m^3]$')
    ax.set_xticks(xticks_pos)
    ax.set_xticklabels(list(data_dict.keys()))

    labels = []
    for q_i in quantiles:
        ind = Indicator.quantile(q_i)
        labels.append(ind.indicator_label)

    ax.legend(handles=err_bars, labels=labels, loc=(0.01, 0.7))

    #ax.set_yscale("log")
    plt.savefig("quantiles.pdf")
    plt.show()
