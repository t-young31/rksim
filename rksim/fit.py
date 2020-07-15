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
    init_concs = [system.network.nodes[i]['c0'] for i in system.network.nodes]

    # Array of times along which the ODE will be solved
    times = np.linspace(0.0, data.get_max_time(), num=1000)

    # Set the rate constants
    # TODO optimise ks
    system.set_init_ks()

    #          dy/dt              y0          t      in notation in the docs
    concs = odeint(system.derivative, init_concs, times)

    # Set the time series for all components in the system
    system.set_simulated(concs, times)
    return None
