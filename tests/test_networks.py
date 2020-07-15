from rksim.reactions import Irreversible
from rksim.species import Reactant, Product
from rksim.systems import System


def test_reaction_network():

    reaction = Irreversible(Reactant(name='R'), Product(name='P'))
    system = System(reaction)

    assert system.network is not None

    # System only has two species and one edge (reaction) between them
    assert system.network.number_of_nodes() == 2
    assert system.network.number_of_edges() == 1


def test_reaction_network_addition():

    reaction = Irreversible(Reactant(name='A'), Reactant(name='B'),
                            Product(name='P'))
    system = System(reaction)

    assert system.network is not None

    # System should have 3 nodes
    assert system.network.number_of_nodes() == 3

    # and 4 edges, one of which is not a reaction but a addition edge
    assert system.network.number_of_edges() == 4
    assert len([edge for edge in system.network.edges
                if system.network.edges[edge]['add'] is True])


def test_set_neighbours():

    reaction = Irreversible(Reactant(name='A'), Reactant(name='B'),
                            Product(name='P'))
    system = System(reaction)

    index_a = next(i for i in system.network.nodes if
                   system.network.nodes[i]['name'] == 'A')

    assert system.network.nodes[index_a]['n_neighbours'] == 1
