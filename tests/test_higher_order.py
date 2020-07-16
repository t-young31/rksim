from rksim.reactions import Irreversible
from rksim.systems import System
from rksim.species import Reactant, Product
from rksim.data import Data, extract_data
from rksim.networks import node_name_to_index
import numpy as np
import os

here = os.path.abspath(os.path.dirname(__file__))


def test_simple():
    # 2R -> P
    reaction = Irreversible(Reactant('R'), Reactant('R'), Product('P'))
    system = System(reaction)

    r_index = node_name_to_index('R', system.network)
    p_index = node_name_to_index('P', system.network)

    # Change in stoichiometry is 2 -> 1 along this reaction edge
    assert system.network.edges[(r_index, p_index)]['sto'] == (2, 1)

    data_path = os.path.join(here, 'higher_order_data', 'second_order.csv')
    data = Data()
    data += extract_data(filename=data_path, names=['R', 'P'])

    data.fit(system)

    # system.plot(name='2nd_order')

    # Should be able to fit the second order data
    assert np.abs(system.get_mse()) < 1E-3

    # Check the initial concentrations and rate constants are as expected.
    # Generated with:  r0 = 0.6,  k = 2
    assert system.network.nodes[r_index]['c0'] == 0.6
    assert np.abs(system.get_rate_constant('R', 'P') - 2.0) < 1E-2
