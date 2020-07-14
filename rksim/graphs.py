import networkx as nx


def make_network(system):
    """
    Make a reaction network for a system of equations

    :param system: (rksim.systems.System)
    """
    graph = nx.Graph()

    # Add all the species to the graph as nodes
    for i, species in enumerate(system.species()):
        graph.add_node(i, name=species.name)

    for reaction in system.reactions:
        pass

    return graph
