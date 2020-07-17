from rksim.species import Reactant, Product
from rksim.reactions import Reversible, Irreversible
from rksim.data import extract_data
from rksim.systems import System
from rksim.data import Data
import numpy as np
import os

here = os.path.abspath(os.path.dirname(__file__))


def test_simple():
    # R -> P

    system = System(Irreversible(Reactant('R'), Product('P')))
    system.set_rate_constants(k=1.0)

    data_path = os.path.join(here, 'simple_data', 'first_order.csv')
    data = Data()

    data += extract_data(filename=data_path, names=['P', 'R'])

    # Fit the system to the data
    data.fit(system, optimise=False)

    # Data is generated with k = 1.0 so the simulated should be essentially
    # exact
    mse = system.mse()
    assert np.abs(mse) < 1E-3

    # -----------------------------------------------------------------------
    # Swapping rate constants should afford a larger error

    system.set_rate_constants(k=2.0)
    data.fit(system, optimise=False)
    # system.plot(name='simple_non_fitted')
    assert np.abs(system.mse()) > 1

    # ----------------------------------------------------------------------
    # Swapping rate constants should a small error when optimised

    system.set_rate_constants(k=2.0)
    data.fit(system, optimise=True)
    # system.plot(name='simple_fitted')

    assert np.abs(system.mse()) < 1E-3


def test_simple_reversible():
    # R <-> P
    system = System(Reversible(Reactant('R'), Product('P')))
    system.set_rate_constants(k=1.0)

    data_path = os.path.join(here, 'simple_data', 'first_order_reversible.csv')
    data = Data()

    # .csv has data in the order time, [P], [R]
    data += extract_data(filename=data_path, names=['P', 'R'])

    # Fit the system to the data, should be able to fit it well
    data.fit(system, optimise=True)
    # system.plot(name='simple_reversible')

    mse = system.mse()
    assert np.abs(mse) < 1E-3

    # Data generated with kf = 0.5,  kb = 0.1
    kf = system.rate_constant('R', 'P')
    assert np.abs(kf - 0.5) < 1E-3

    kb = system.rate_constant('P', 'R')
    assert np.abs(kb - 0.1) < 1E-3


def test_simple_alt():
    # B <-- A --> C
    system = System(Irreversible(Reactant('A'), Product('B')),
                    Irreversible(Reactant('A'), Product('C')))
    system.set_rate_constants(k=1.0)

    data_path = os.path.join(here, 'simple_data', 'double_first_order.csv')
    data = Data()
    data += extract_data(filename=data_path, names=['A', 'B', 'C'])

    data.fit(system)
    # system.plot(name='double_first_order')

    # Should be able to fit these data well
    mse = system.mse()
    assert np.abs(mse) < 1E-3

    # Data generated with: k1 = 0.4, k2 = 0.8
    assert np.abs(system.rate_constant('A', 'B') - 0.4) < 1E-2
    assert np.abs(system.rate_constant('A', 'C') - 0.8) < 1E-2
