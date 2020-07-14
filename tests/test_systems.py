from rksim.systems import *
from rksim.data import Data, extract_data, TimeSeries
import os

here = os.path.abspath(os.path.dirname(__file__))


def test_simple_system():

    system = System(Irreversible(Reactant(name='R'),
                                 Product(name='P')))

    data_path = os.path.join(here, 'data', 'first_order.csv')

    data = Data()
    data += extract_data(data_path, names=['R', 'P'])

    data.assign_series(system=system)

    reaction = system.reactions[0]
    reactant = reaction.reactants[0]

    assert reactant.name == 'R'
    assert reactant.time_series is not None
    assert isinstance(reactant.time_series, TimeSeries)

    product = reaction.products[0]

    assert product.name == 'P'
    assert product.time_series is not None
    assert isinstance(product.time_series, TimeSeries)


def test_reaction():

    reaction = Reaction(Reactant(name='R'), Reactant(name='R'),
                        Product(name='P'))

    # Reaction on;y has two compoents, one js just repeated
    assert len(reaction.components) == 2

    # R^2 -> P
    assert reaction.components[0].order == 2
    assert reaction.components[1].order == 1
