from rksim.species import Reactant, Product
from rksim.reactions import Irreversible
from rksim.data import extract_data
from rksim.systems import System
from rksim.data import Data
import numpy as np
import os

here = os.path.abspath(os.path.dirname(__file__))


def test_simple():

    system = System(Irreversible(Reactant('R'), Product('P')))
    system.set_rate_constants(k=1.0)

    data_path = os.path.join(here, 'data', 'first_order.csv')
    data = Data()

    data += extract_data(filename=data_path, names=['P', 'R'])

    # Fit the system to the data
    data.fit(system, optimise=False)

    # Data is generated with k = 1.0 so the simulated should be essentially
    # exact
    mse = system.get_mse()
    assert np.abs(mse) < 1E-3

    # ----------------------------------------------------
    # Swapping rate constants should afford a larger error

    system.set_rate_constants(k=2.0)
    data.fit(system, optimise=False)

    assert np.abs(system.get_mse()) > 1
