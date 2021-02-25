from typing import Union
from scipy.integrate import odeint
from scipy.optimize import minimize
import numpy as np


def fit(data, system, optimise: Union[list, bool], max_time=None):
    """
    Fit a system of equations to self.series. If not specified
    (system=None) then all reasonable fits will be attempted. This
    function will append to self.fits

    :param data: (rksim.data.Data)

    :param optimise: (bool | list) Whether to optimise the rate constant(s), if
                     a list then will optimise the indexes of the rate constants
                     defined in the list

    :param system: (rksim.system.System) equation system that these
                   data will be fit to
    """
    system.set_stoichiometries()

    init_concs = [system.network.nodes[i]['c0'] for i in system.network.nodes]

    # Array of times along which the ODE will be solved
    if max_time is None and data is not None:
        max_time = data.max_time()

    times = np.linspace(0.0, max_time, num=10000)

    if optimise is not False:
        # Minimise the difference between the simulated and observed data wrt.
        # the rate constants and set the optimised values

        k_idxs_to_opt = [idx for idx in range(len(system))]
        init_ks = system.rate_constants()

        if type(optimise) is list:
            k_idxs_to_opt = np.array([idx for idx in optimise], dtype=int)
            init_ks = init_ks[k_idxs_to_opt]

        result = minimize(mse,
                          x0=init_ks,
                          method='BFGS',
                          args=(system, times, init_concs, k_idxs_to_opt))

        # Rate constants must all be positive (abs(k)) is minimised in mse()
        for i, idx in enumerate(k_idxs_to_opt):
            system.set_rate_constant(idx, k=np.abs(result.x)[i])

        return None

    #                    dy/dt            y0        t    in scipy doc notation
    concs = odeint(system.derivative, init_concs, times)

    # Set the time series for all components in the system
    system.set_simulated(concs, times)
    return None


def mse(rate_constants, system, times, init_concs, k_idxs_to_opt):
    """Calculate the rel mean squared error for a system with a set of ks"""

    # Set only the rate constants that are going to be optimised with the
    # new values
    for i, idx in enumerate(k_idxs_to_opt):
        system.set_rate_constant(idx, k=rate_constants[i])

    # Integrate the time forward and set the time series
    concs = odeint(system.derivative, init_concs, times)
    system.set_simulated(concs, times)

    return system.mse(relative=True)
