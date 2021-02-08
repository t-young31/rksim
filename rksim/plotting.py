import matplotlib.pyplot as plt
from rksim.exceptions import PlottingFailed


def plot(generator, name=None, dpi=400):
    """Plot a series of times series with matplotlib"""
    # Plot with consistent colors
    cm = plt.get_cmap('tab10')

    for i, item in enumerate(generator):
        color = cm.colors[i % 10]

        # Plot a time series directly
        if hasattr(item, 'times'):
            plt.scatter(item.times, item.concentrations,  label=f'{item.name}',
                        marker='o', s=5, color=color)
            continue

        # Must otherwise be a series
        if not hasattr(item, 'series'):
            raise PlottingFailed('Can only plot TimeSeries or Species')

        # Plot a time series (expt or simulated from a rksim.species.Species)
        if item.series is not None:
            plt.scatter(item.series.times, item.series.concentrations,
                        label=f'{item.series.name}',
                        marker='o', s=5, color=color)

        if item.simulated_series is not None:
            series = item.simulated_series

            plt.plot(series.times, series.concentrations,
                     label=f'{series.name} simulated',
                     lw=1.0, ls='--', color=color)

    # Legend and axis labels
    plt.legend()
    plt.ylabel('Concentration / mol dm$^{-3}$')
    plt.xlabel('Time / s')

    return show_or_plot(name, dpi)


def show_or_plot(name=None, dpi=400):
    """Show the plt if name=None otherwise save it"""
    plt.tight_layout()

    if name is not None:
        plt.savefig(f'{name}.png', dpi=dpi)

    else:
        plt.show()

    plt.close()
    return None
