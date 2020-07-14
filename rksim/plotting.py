import matplotlib.pyplot as plt


def plot(expt=None, sim=None, name=None, dpi=400):
    """Plot a series of times series with matplotlib"""

    # Plot any experimental points
    if expt is not None:
        for time_series in expt:
            plt.scatter(time_series.times, time_series.concentrations,
                        label=f'{time_series.name}', marker='o', s=3)

    # Plot any simulated data
    if sim is not None:
        for time_series in sim:
            plt.plot(time_series.times, time_series.concentrations,
                     label=f'{time_series.name} simulated', lw=1.0, ls='--')

    # Legend and axis labels
    plt.legend()
    plt.ylabel('Concentration / mol dm$^{-3}$')
    plt.xlabel('Time / s')

    if name is not None:
        plt.savefig(name, dpi=dpi)

    else:
        plt.show()

    plt.close()
    return None
