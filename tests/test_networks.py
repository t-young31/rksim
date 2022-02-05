from rksim.reactions import ReversibleReaction, IrreversibleReaction
from rksim.species import Reactant, Product
from rksim.systems import System
import os


def test_reaction_network():

    reaction = IrreversibleReaction(Reactant(name='R'), Product(name='P'))
    system = System(reaction)

    assert system.network is not None

    # System only has two species
    assert system.network.number_of_nodes() == 2


def test_plotting_network():
    reaction = IrreversibleReaction(Reactant(name='R'), Product(name='P'))
    system = System(reaction)

    system.network.plot(name='test_network')
    assert os.path.exists('test_network.pdf')
    os.remove('test_network.pdf')


def test_reaction_network_addition():

    reaction = IrreversibleReaction(Reactant(name='A'), Reactant(name='B'),
                                    Product(name='P'))
    system = System(reaction)

    assert system.network is not None

    # System should have 3 nodes
    assert system.network.number_of_nodes() == 3


def test_rate_constants():
    # A + B -> P
    reaction = IrreversibleReaction(Reactant('A'), Reactant('B'), Product('P'))
    system = System(reaction)

    system.set_rate_constants(ks=[0.5])

    ks = system.rate_constants()
    assert all(k == 0.5 for k in ks)

    # A + B -> C + D
    reaction = IrreversibleReaction(Reactant('A'), Reactant('B'),
                                    Product('C'), Product('D'))
    system = System(reaction)
    print(system)


def test_reversible_rate_constants():
    # A + B <-> P
    reaction = ReversibleReaction(Reactant('A'), Reactant('B'), Product('P'))
    system = System(reaction)

    # Should have forward and backward rate constants
    assert len(system.rate_constants()) == 2

    # A + B <-> P + D
    system = System(ReversibleReaction(Reactant('A'), Reactant('B'),
                                       Product('P'), Product('D')))
    assert len(system.rate_constants()) == 2

    # A + B <-> C <-> D
    system = System(ReversibleReaction(Reactant('A'), Product('C')),
                    ReversibleReaction(Reactant('C'), Product('D')))
    assert len(system.rate_constants()) == 4

    # A + B -> C <-> D + E -> F
    system = System(IrreversibleReaction(Reactant('A'), Reactant('B'), Product('C')),
                    ReversibleReaction(Reactant('C'), Product('D'), Product('E')),
                    IrreversibleReaction(Reactant('E'), Product('F')))

    assert len(system.rate_constants()) == 4


def test_loop():
    # E + S <--> ES -> P + E
    # Classic Michaelis-Menten kinetics

    system = System(ReversibleReaction(Reactant('S'), Reactant('E'), Product('ES')),
                    IrreversibleReaction(Reactant('ES'), Product('E'), Product('P')))

    # Should be 3 distinct rate constants for this system
    assert len(system.rate_constants()) == 3
