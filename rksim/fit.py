from scipy.integrate import odeint
from scipy.optimize import minimize
import numpy as np


def fit(data, system, optimise):
    """
    Fit a system of equations to self.series. If not specified
    (system=None) then all reasonable fits will be attempted. This
    function will append to self.fits

    :param data: (rksim.data.Data)

    :param optimise: (bool) Whether to optimise the rate constant(s)

    :param system: (rksim.system.System) equation system that these
                   data will be fit to
    """
    init_concs = [system.network.nodes[i]['c0'] for i in system.network.nodes]

    # Array of times along which the ODE will be solved
    times = np.linspace(0.0, data.max_time(), num=1000)

    if optimise:
        optimise_rate_constants(system, times, init_concs)

    #                    dy/dt            y0        t    in scipy doc notation
    concs = odeint(system.derivative, init_concs, times)

    # Set the time series for all components in the system
    system.set_simulated(concs, times)
    return None


def optimise_rate_constants(system, times, init_concs):
    """
    Optimise the rate constants in a system of equations to minimise
    the mean squared error between the simulated and experimental data

    :param system: (rksim.system.System) equation system that these
                   data will be fit to

    :param times: (np.ndarray) Times (s) over which to solve the ODEs

    :param init_concs: (np.ndarray) Initial concentrations (mol dm^-3)
    """

    init_ks = system.rate_constants()

    # Minimise the error and set the optimise rate constants
    result = minimize(mse, x0=init_ks, args=(system, times, init_concs))
    system.set_rate_constants(result.x)

    return None


def mse(rate_constants, system, times, init_concs):
    """Calculate the mean squared error for a system with a set of ks"""

    system.set_rate_constants(rate_constants)

    # Integrate the time forward and set the time series
    concs = odeint(system.derivative, init_concs, times)
    system.set_simulated(concs, times)

    return system.mse()
