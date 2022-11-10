import numpy as np
import matplotlib.pyplot as plt

def __calculate_plot_data(variogram):
    # get the parameters
    _bins = variogram.bins
    _exp = variogram.experimental
    x = np.linspace(0, np.nanmax(_bins), 100)

    # apply the model
    y = variogram.transform(x)

    # handle the relative experimental variogram
    if variogram.normalized:
        _bins /= np.nanmax(_bins)
        y /= np.max(_exp)
        _exp /= np.nanmax(_exp)
        x /= np.nanmax(x)

    return x, y, _bins, _exp


def matplotlib_variogram_plot(
    variograms,
    axes=None,
    grid=True,
    show=True,
    hist=True,
    title=None
):
    # do the plotting
    if axes is None:
        if hist:
            fig = plt.figure(figsize=(8, 5))
            ax1 = plt.subplot2grid((5, 1), (1, 0), rowspan=4)
            ax2 = plt.subplot2grid((5, 1), (0, 0), sharex=ax1)
            fig.subplots_adjust(hspace=0)
        else:
            fig, ax1 = plt.subplots(1, 1, figsize=(8, 4))
            ax2 = None
    elif isinstance(axes, (list, tuple, np.ndarray)):
        ax1, ax2 = axes
        fig = ax1.get_figure()
    else:
        ax1 = axes
        ax2 = None
        fig = ax1.get_figure()

    colors = ['blue', 'green', 'red']


    for index, (variogram, label, subdir) in enumerate(variograms):
        # get the plotting data
        x, y, _bins, _exp = __calculate_plot_data(variogram)

        print("x ", x)
        print("y ", y)

        # ------------------------
        # plot Variograms model
        ax1.plot(_bins, _exp, ".", color=colors[index])
        ax1.plot(x, y, '-', color=colors[index], label=label)

        # ax limits
        if variogram.normalized:
            ax1.set_xlim([0, 1.05])
            ax1.set_ylim([0, 1.05])

        # grid settings
        if grid:
            ax1.grid(False)
            ax1.vlines(
                _bins,
                *ax1.axes.get_ybound(),
                colors=(.85, .85, .85),
                linestyles='dashed'
            )

        # always print error bars above grid
        conf = variogram._experimental_conf_interval
        if conf is not None:
            lo = conf[:, 1] - conf[:, 0]
            up = conf[:, 2] - conf[:, 1]
            yerr = np.column_stack((lo, up)).T
            ax1.errorbar(_bins, _exp, fmt='.b', yerr=yerr)

        # annotation
        ax1.axes.set_ylabel('semivariance (%s)' % variogram._estimator.__name__)
        ax1.axes.set_xlabel('Lag (-)')

        # ------------------------
        # plot histogram
        if ax2 is not None and hist:
            # calc the histogram
            _count = np.fromiter(
                (g.size for g in variogram.lag_classes()), dtype=int
            )

            # set the sum of hist bar widths to 70% of the x-axis space
            w = (np.max(_bins) * 0.7) / len(_count)

            # plot
            ax2.bar(_bins, _count, width=w, align='center', color='red')

            # adjust
            plt.setp(ax2.axes.get_xticklabels(), visible=False)
            ax2.axes.set_yticks(ax2.axes.get_yticks()[1:])

            # need a grid?
            if grid:  # pragma: no cover
                ax2.grid(False)
                ax2.vlines(
                    _bins,
                    *ax2.axes.get_ybound(),
                    colors=(.85, .85, .85),
                    linestyles='dashed'
                )

            # anotate
            ax2.axes.set_ylabel('N')

    fig.legend()

    # show the figure
    if show:  # pragma: no cover
        fig.show()

    fig.savefig(title)

    return fig