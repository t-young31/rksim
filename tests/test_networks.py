from rksim.systems import *


def test_reaction_network():

    reaction = Irreversible(Reactant(name='R'), Product(name='P'))
    system = System(reaction)

    assert system.network is not None

    # System only has two species and one edge (reaction) between them
    assert system.network.number_of_nodes() == 2
    assert system.network.number_of_edges() == 1
