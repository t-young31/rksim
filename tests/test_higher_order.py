from rksim.reactions import Irreversible
from rksim.systems import System
from rksim.species import Reactant, Product
from rksim.data import Data, extract_data
import numpy as np
import os

here = os.path.abspath(os.path.dirname(__file__))


def test_simple():
    # 2R -> P
    reaction = Irreversible(Reactant('R'), Reactant('R'), Product('P'))
    system = System(reaction)

    r_index = system.network.node_mapping['R']

    data_path = os.path.join(here, 'higher_order_data', 'second_order.csv')
    data = Data()
    data += extract_data(filename=data_path, names=['R', 'P'])

    data.fit(system)
    # system.plot(name='2nd_order')

    # Should be able to fit the second order data
    assert np.abs(system.mse()) < 1E-3

    # Check the initial concentrations and rate constants are as expected.
    # Generated with:  r0 = 0.6,  k = 2
    assert system.network.nodes[r_index]['c0'] == 0.6
    assert np.abs(system.rate_constant('R', 'P') - 2.0) < 1E-2


def test_ab_reaction():

    reaction = Irreversible(Reactant('A'), Reactant('B'),
                            Product('P'))

    system = System(reaction)

    data_path = os.path.join(here, 'higher_order_data', 'ab.csv')
    data = Data()
    data += extract_data(filename=data_path, names=['A', 'B', 'P'])

    data.fit(system)
    # system.plot(name='a_b')

    # Should be able to fit the data well
    assert np.abs(system.mse()) < 1E-2

    # Data generated with k = 1.7 s^-2 mol^-1 dm^3
    assert np.abs(system.rate_constant('A', 'B',  'P') - 1.7) < 1E-2

    ks = system.rate_constants()
    assert len(ks) == 1
    assert np.abs(ks[0] - 1.7) < 1E-2
