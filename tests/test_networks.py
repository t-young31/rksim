from rksim.reactions import Reversible, Irreversible
from rksim.species import Reactant, Product
from rksim.systems import System
from rksim.networks import node_name_to_index


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


def test_edge_mapping():
    # A + B -> P
    reaction = Irreversible(Reactant('A'), Reactant('B'), Product('P'))
    system = System(reaction)

    node_a = node_name_to_index('A', graph=system.network)
    node_b = node_name_to_index('B', graph=system.network)
    node_p = node_name_to_index('P', graph=system.network)

    key = next(iter(system.network.edge_mapping.keys()))

    assert (node_a, node_p) in system.network.edge_mapping[key]     # A-P
    assert (node_b, node_p) in system.network.edge_mapping[key]     # B-P

    assert (node_a, node_b) not in system.network.edge_mapping[key]  # A-B
    assert (node_p, node_p) not in system.network.edge_mapping[key]  # P-P


def test_rate_constants():
    # A + B -> P
    reaction = Irreversible(Reactant('A'), Reactant('B'), Product('P'))
    system = System(reaction)

    system.set_rate_constants(ks=[0.5])

    for (i, j) in system.network.edges:
        if system.network[i][j]['add'] is False:
            assert system.network[i][j]['k'] == 0.5

    ks = system.get_rate_constants()
    assert all(k == 0.5 for k in ks)

    # A + B -> C + D
    reaction = Irreversible(Reactant('A'), Reactant('B'),
                            Product('C'), Product('D'))
    system = System(reaction)
    print(system)

    print(system.network.edge_mapping)
    assert len(system.get_rate_constants()) == 1


def test_reversible_rate_constants():
    # A + B <-> P
    reaction = Reversible(Reactant('A'), Reactant('B'), Product('P'))
    system = System(reaction)

    # Should have forward and backward rate constants
    assert len(system.get_rate_constants()) == 2

    # A + B <-> P + D
    system = System(Reversible(Reactant('A'), Reactant('B'),
                               Product('P'), Product('D')))
    assert len(system.get_rate_constants()) == 2

    # A + B <-> C <-> D
    system = System(Reversible(Reactant('A'), Product('C')),
                    Reversible(Reactant('C'), Product('D')))
    assert len(system.get_rate_constants()) == 4

    # A + B -> C <-> D + E -> F
    system = System(Irreversible(Reactant('A'), Reactant('B'), Product('C')),
                    Reversible(Reactant('C'), Product('D'), Product('E')),
                    Irreversible(Reactant('E'), Product('F')))

    assert len(system.get_rate_constants()) == 4
