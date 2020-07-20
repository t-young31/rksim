from rksim.systems import System
from rksim.reactions import Irreversible, Reversible, Reaction
from rksim.species import Reactant, Product
from rksim.exceptions import CannotSetAttribute, CannotGetAttribute
from rksim.data import Data, extract_data, TimeSeries
import pytest
import numpy as np
import os

here = os.path.abspath(os.path.dirname(__file__))


def test_simple_system():

    system = System(Irreversible(Reactant(name='R'),
                                 Product(name='P')))

    data_path = os.path.join(here, 'simple_data', 'first_order.csv')

    data = Data()
    data += extract_data(data_path, names=['R', 'P'])

    data.assign(system=system)

    for species in system.species():
        assert species.name in ['R', 'P']

        assert species.series is not None
        assert isinstance(species.series, TimeSeries)


def test_mse():
    system = System(Irreversible(Reactant('R'), Product('P')))

    # System not assigned to any data shouldn't have any mean squared error
    assert system.mse() == 0


def test_reaction():

    reaction = Reaction(Reactant(name='R'), Reactant(name='R'),
                        Product(name='P'))

    # Reaction on;y has two components, R is repeated
    assert len(reaction.components) == 2

    # 2R -> P
    for species in reaction.components:
        if species.name == 'R':
            assert species.stoichiometry == 2

        if species.name == 'P':
            assert species.stoichiometry == 1


def test_species():
    """
    A system can be initialised from reactions with the same species
    but they should not be present more than once in system.species()
    """
    reaction1 = Irreversible(Reactant(name='R'), Product(name='P'))
    reaction2 = Irreversible(Reactant(name='A'), Product(name='P'))

    system = System(reaction1, reaction2)

    names = []

    for species in system.species():
        names.append(species.name)

    assert len(names) == 3          # Should only have 3 species
    assert 'R' in names
    assert 'A' in names
    assert 'P' in names


def test_setting_k():
    system = System(Irreversible(Reactant(name='R'),
                                 Product(name='P')))
    system.set_rate_constants(k=10.0)

    for edge in system.network.edges:
        assert system.network.edges[edge]['k'] == 10


def test_fit():
    system = System(Irreversible(Reactant('R'), Product('P')))
    system.set_rate_constants(k=1.2)

    data_path = os.path.join(here, 'simple_data', 'first_order.csv')

    data = Data()
    data += extract_data(data_path, names=['P', 'R'])

    # Can do system.fit(data) or data.fit(system)
    system.fit(data)

    # Should be able to fit these data
    assert np.abs(system.mse()) < 1E-2

    # Data generated with k=1.0 s^-1
    assert np.abs(system.rate_constant('R', 'P') - 1.0) < 1E-2


def test_derivative1():
    # R -> P
    system = System(Irreversible(Reactant(name='R'),
                                 Product(name='P')))
    # Set some reasonable rates
    system.set_rate_constants(k=1.0)

    # Concentration of P, R as components in a reaction are
    # sorted alphabetically
    concs = [0.0, 1.0]

    dpdt, drdt = system.derivative(concs)

    expected_dpdt = 1.0
    assert np.abs(dpdt - expected_dpdt) < 1E-8

    expected_drdt = -1.0
    assert np.abs(drdt - expected_drdt) < 1E-8


def test_derivative2():
    # A + B -> C
    system = System(Irreversible(Reactant(name='A'), Reactant(name='B'),
                                 Product(name='C')))
    k = 1.0
    system.set_rate_constants(k=k)

    possible_concentrations = [[2.0, 1.0, 0.0],
                               [0.0, 1.0, 0.0],
                               [10.0, 0.0, 1.0],
                               [1.0, 1.0, 0.0]]

    for c_a, c_b, c_c in possible_concentrations:

        dadt, dbdt, dcdt = system.derivative([c_a, c_b, c_c])

        assert np.abs(dadt - (-k * c_a * c_b)) < 1E-8
        assert np.abs(dbdt - (-k * c_a * c_b)) < 1E-8
        assert np.abs(dcdt - (+k * c_a * c_b)) < 1E-8


def test_derivative3():
    # A -> B + C
    system = System(Irreversible(Reactant(name='A'),
                                 Product(name='B'), Product(name='C')))
    k = 1.0
    system.set_rate_constants(k=k)
    c_a, c_b, c_c = [2.0, 1.0, 0.0]

    dadt, dbdt, dcdt = system.derivative([c_a, c_b, c_c])

    assert np.abs(dadt - (-k * c_a)) < 1E-8
    assert np.abs(dbdt - (+k * c_a)) < 1E-8
    assert np.abs(dcdt - (+k * c_a)) < 1E-8


def test_derivative4():
    # R  -> P    (reversible)
    #    <-
    system = System(Reversible(Reactant(name='R'),
                               Product(name='P')))

    kf = kb = 1.0
    system.set_rate_constants(k=kf)

    #                          [P]   [R]
    possible_concentrations = [[1.0, 0.0],
                               [0.5, 0.5],
                               [0.1, 0.9]]

    for c_p, c_r in possible_concentrations:

        dpdt, drdt = system.derivative([c_p, c_r])

        assert np.abs(dpdt - (kf * c_r - kb * c_p)) < 1E-8
        assert np.abs(drdt - (-kf * c_r + kb * c_p)) < 1E-8


def test_derivative5():
    # A + B -> C + D
    #       <-
    system = System(Reversible(Reactant(name='A'), Reactant(name='B'),
                               Product(name='C'), Product(name='D')))
    k = 1.0
    system.set_rate_constants(k=k)
    c_a, c_b, c_c, c_d = [2.0, 1.0, 0.2, 0.1]

    dadt, dbdt, dcdt, dddt = system.derivative([c_a, c_b, c_c, c_d])

    assert np.abs(dadt - (-k * c_a * c_b + k * c_c * c_d)) < 1E-8
    assert np.abs(dbdt - (-k * c_a * c_b + k * c_c * c_d)) < 1E-8
    assert np.abs(dcdt - (+k * c_a * c_b - k * c_c * c_d)) < 1E-8
    assert np.abs(dddt - (+k * c_a * c_b - k * c_c * c_d)) < 1E-8


def test_derivative6():
    # A
    #  \
    #   C -> D
    #  /
    # B
    system = System(Irreversible(Reactant('A'), Product('C')),
                    Irreversible(Reactant('B'), Product('C')),
                    Irreversible(Reactant('C'), Product('D')))
    k = 1.0
    system.set_rate_constants(k=k)

    #                     A    C   B    D
    c_a, c_c, c_b, c_d = [1.0, 2.0, 0.5, 0.0]

    dadt, dcdt, dbdt, dddt = system.derivative([c_a, c_c, c_b, c_d])

    assert np.abs(dadt - (-k * c_a)) < 1E-8
    assert np.abs(dbdt - (-k * c_b)) < 1E-8
    assert np.abs(dcdt - (k * c_a + k * c_b - k * c_c)) < 1E-8
    assert np.abs(dddt - (k * c_c)) < 1E-8


def test_derivative7():
    # 2A -> B
    system = System(Irreversible(Reactant('A'), Reactant('A'),
                                 Product('B')))
    k = 1.0
    system.set_rate_constants(k=k)

    # Check the graph has the correct order in this component
    a_index = system.network.node_mapping['A']
    b_index = system.network.node_mapping['B']
    assert system.network.edges[(a_index, b_index)]['sto'] == (2, 1)

    c_a, c_b = [2.0, 0.5]

    dadt, dbdt = system.derivative([c_a, c_b])

    assert np.abs(dadt - (-2 * k * c_a**2)) < 1E-8
    assert np.abs(dbdt - (k * c_a**2)) < 1E-8


def test_derivative8():
    # A -> 2B
    system = System(Irreversible(Reactant('A'),
                                 Product('B'), Product('B')))
    k = 1.0
    system.set_rate_constants(k=k)

    c_a, c_b = [2.0, 0.5]

    dadt, dbdt = system.derivative([c_a, c_b])

    assert np.abs(dadt - (-k * c_a)) < 1E-8
    assert np.abs(dbdt - (+2.0 * k * c_a)) < 1E-8


def test_derivative9():
    # 2A -> B, A -> C

    reaction1 = Irreversible(Reactant('A'), Reactant('A'), Product('B'))
    reaction2 = Irreversible(Reactant('A'), Product('C'))

    system = System(reaction1, reaction2)

    k1, k2 = 0.3, 0.5
    system.set_rate_constants([k1, k2])

    c_a, c_b, c_c = 0.7, 0.0, 0.0

    dadt, dbdt, dcdt = system.derivative([c_a, c_b, c_c])

    # expecting: d[A]/dt = -2k1[A]^2 - k2[A]
    assert np.abs(dadt - (-2 * k1 * c_a**2 - k2 * c_a)) < 1E-8


def test_set_init_conc():
    # R -> P
    system = System(Irreversible(Reactant('R'), Product('P')))

    system.set_initial_concentration(name='R', c=2.0)

    for node in system.network.nodes:

        if system.network.nodes[node]['name'] == 'R':
            assert system.network.nodes[node]['c0'] == 2.0

    # Cannot set the concentration for a species that doesn't exist
    with pytest.raises(CannotSetAttribute):
        system.set_initial_concentration(name='A', c=2.0)


def test_set_rate_constant():
    # R -> P
    system = System(Irreversible(Reactant('R'), Product('P')))

    system.set_rate_constant('R', 'P', k=2.0)

    for edge in system.network.edges:
        assert system.network.edges[edge]['k'] == 2.0

    # No reverse reaction
    with pytest.raises(CannotSetAttribute):
        system.set_rate_constant('P', 'R', k=2.0)


def test_get_rate_constant():

    system = System(Irreversible(Reactant('R'), Product('P')))

    assert system.rate_constant('R', 'P') is not None

    with pytest.raises(CannotGetAttribute):
        system.rate_constant('X', 'Y')
