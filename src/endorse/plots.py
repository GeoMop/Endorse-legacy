from typing import *
import matplotlib.pyplot as plt
from matplotlib import ticker
from .indicator import IndicatorFn

import numpy as np
import pandas as pd


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
    # how to use promille
    # plt.rcParams['text.latex.preamble'] = [r"\usepackage{wasysym}"]
    # \textperthousand
    plt.rcParams.update({'font.size': 12})

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

