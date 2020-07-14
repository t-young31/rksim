from scipy.integrate import odeint
import numpy as np


def fit(data, system):
    """
    Fit a system of equations to self.series. If not specified
    (system=None) then all reasonable fits will be attempted. This
    function will append to self.fits

    :param data: (rksim.data.Data)

    :param system: (rksim.system.System) equation system that these
                   data will be fit to
    """
    init_concs = get_initial_concentrations(system)

    # Array of times along which the ODE will be solved
    times = np.linspace(0.0, data.get_max_time(), num=1000)

    # Set the rate constants
    system.set_init_ks()

    #          dy/dt              y0          t      in notation in the docs
    concs = odeint(system.derivative, init_concs, times)

    # Set the time series for all components in the system
    system.set_simulated(concs, times)
    return None


def get_initial_concentrations(system):
    """Get the initial concentrations for all components in a system

    If there is not a specified time series for a component then the initial
    concentration will be set to zero.

    :param system: (rksim.system.System)
    """

    init_concs = []

    for species in system.species():
        if species.time_series is not None:
            # Take the first concentration from the time series
            c0 = species.time_series.concentrations[0]

        else:
            c0 = 0.0

        if c0 == 0:       # Bias positive propagation
            c0 += 1E-8

        init_concs.append(c0)

    return init_concs
