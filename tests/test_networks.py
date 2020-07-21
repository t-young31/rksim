from rksim.reactions import Reversible, Irreversible
from rksim.species import Reactant, Product
from rksim.systems import System
import os


def test_reaction_network():

    reaction = Irreversible(Reactant(name='R'), Product(name='P'))
    system = System(reaction)

    assert system.network is not None

    # System only has two species
    assert system.network.number_of_nodes() == 2


def test_plotting_network():
    reaction = Irreversible(Reactant(name='R'), Product(name='P'))
    system = System(reaction)

    system.network.plot(name='test_network')
    assert os.path.exists('test_network.png')
    os.remove('test_network.png')


def test_reaction_network_addition():

    reaction = Irreversible(Reactant(name='A'), Reactant(name='B'),
                            Product(name='P'))
    system = System(reaction)

    assert system.network is not None

    # System should have 3 nodes
    assert system.network.number_of_nodes() == 3


def test_rate_constants():
    # A + B -> P
    reaction = Irreversible(Reactant('A'), Reactant('B'), Product('P'))
    system = System(reaction)

    system.set_rate_constants(ks=[0.5])

    ks = system.rate_constants()
    assert all(k == 0.5 for k in ks)

    # A + B -> C + D
    reaction = Irreversible(Reactant('A'), Reactant('B'),
                            Product('C'), Product('D'))
    system = System(reaction)
    print(system)


def test_reversible_rate_constants():
    # A + B <-> P
    reaction = Reversible(Reactant('A'), Reactant('B'), Product('P'))
    system = System(reaction)

    # Should have forward and backward rate constants
    assert len(system.rate_constants()) == 2

    # A + B <-> P + D
    system = System(Reversible(Reactant('A'), Reactant('B'),
                               Product('P'), Product('D')))
    assert len(system.rate_constants()) == 2

    # A + B <-> C <-> D
    system = System(Reversible(Reactant('A'), Product('C')),
                    Reversible(Reactant('C'), Product('D')))
    assert len(system.rate_constants()) == 4

    # A + B -> C <-> D + E -> F
    system = System(Irreversible(Reactant('A'), Reactant('B'), Product('C')),
                    Reversible(Reactant('C'), Product('D'), Product('E')),
                    Irreversible(Reactant('E'), Product('F')))

    assert len(system.rate_constants()) == 4


def test_loop():
    # E + S <--> ES -> P + E
    # Classic Michaelis-Menten kinetics

    system = System(Reversible(Reactant('S'), Reactant('E'), Product('ES')),
                    Irreversible(Reactant('ES'), Product('E'), Product('P')))

    # Should be 3 distinct rate constants for this system
    assert len(system.rate_constants()) == 3
