import matplotlib as mpl
import matplotlib.pyplot as plt
from rksim.exceptions import PlottingFailed

mpl.rcParams['axes.labelsize'] = 15
mpl.rcParams['lines.linewidth'] = 1.5
mpl.rcParams['lines.markersize'] = 5
mpl.rcParams['xtick.labelsize'] = 13
mpl.rcParams['ytick.labelsize'] = 13
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.top'] = True
mpl.rcParams['ytick.right'] = True
mpl.rcParams['axes.linewidth'] = 1.2
plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))


def plot(generator, name=None, time_units='s', conc_units='mol dm$^{-3}$'):
    """Plot a series of times series with matplotlib"""
    # Plot with consistent colors
    cm = plt.get_cmap('tab10')

    for i, item in enumerate(generator):
        color = cm.colors[i % 10]

        # Plot a time series directly
        if hasattr(item, 'times'):
            plt.scatter(item.times, item.concentrations,  label=f'{item.name}',
                        marker='o', s=15, color=color)
            continue

        # Must otherwise be a series
        if not hasattr(item, 'series'):
            raise PlottingFailed('Can only plot TimeSeries or Species')

        # Plot a time series (expt or simulated from a rksim.species.Species)
        if item.series is not None:
            plt.scatter(item.series.times, item.series.concentrations,
                        label=f'{item.series.name}',
                        marker='o', s=15, color=color)

        if item.simulated_series is not None:
            series = item.simulated_series

            plt.plot(series.times, series.concentrations,
                     label=f'{series.name} simulated',
                      ls='--', color=color)

    # Legend and axis labels
    plt.legend()
    conc_units = "" if conc_units is None else f'/ {conc_units}'
    plt.ylabel(f'Concentration {conc_units}')

    time_units = "" if time_units is None else f'/ {time_units}'
    plt.xlabel(f'Time {time_units}')

    return show_or_plot(name)


def show_or_plot(name=None):
    """Show the plt if name=None otherwise save it"""
    plt.tight_layout()

    if name is not None:
        plt.savefig(f'{name}.pdf')

    else:
        plt.show()

    plt.close()
    return None
