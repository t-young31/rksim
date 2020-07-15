import networkx as nx


def make_network(reactions):
    """
    Make a reaction network for a system of reactions


    e.g. For an irreversible A + B -> P reaction

           P
          ^ ^
         /  \
       /     \
      A ----- B
          +

    where the arrowed edges have an associated rate constant attribute (k)
    and the addition edge between A and B has add=True

    :param reactions: (list(rksim.systems.Reaction))
    :return: (networkx.Graph)
    """
    graph = nx.DiGraph()

    # Add all species as nodes to the graph/network
    add_nodes(graph, reactions=reactions)

    for reaction in reactions:
        # Add any + edges between reactants, then products
        add_addition_edges(graph, species_generator=reaction.reactants())
        add_addition_edges(graph, species_generator=reaction.products())

        # Add the reaction arrow with associated rate constant (k) values
        add_reaction_edges(graph, reaction=reaction)

    populate_neighbour_number(graph)

    return graph


def inflow_node(graph, node_index):
    """
    Get a nodes that flows into a node_index

    :param graph: (networkx.Graph)
    :param node_index: (int)
    """
    for j in graph.predecessors(node_index):
        if graph[j][node_index]['add'] is False:
            yield j

    return StopIteration


def outflow_node(graph, node_index):
    """
    Get a nodes that flows from a node_index

    :param graph: (networkx.Graph)
    :param node_index: (int)
    """
    for j in graph.successors(node_index):
        if graph[node_index][j]['add'] is False:
            yield j

    return StopIteration


def populate_neighbour_number(graph):
    """
    Populate the number of neighbours over add edges for a reaction graph

           P
          ^ ^
         /  \
       /     \
      A ----- B
          +

    P_n_neighbours = 0, A_n_neighbours = 1, B_n_neighbours = 1

    :param graph: (networkx.Graph)
    """
    neighbours_dict = {}

    for node in graph.nodes:
        n_neighbours = 0

        for neighbour in graph.neighbors(node):
            if graph[node][neighbour]['add']:
                n_neighbours += 1

        neighbours_dict[node] = n_neighbours

    # Set the number of neighbours for each species in the network/graph
    nx.set_node_attributes(graph, neighbours_dict, name='n_neighbours')
    return None


def neighbours(graph, node_index):
    """
    Get the neighbours along addition edges next to this atom. For example
    a graph

        +
    A ----- B

    where A is indexed 0 and B 1 then get_add_neighbours(graph, 0) -> [1]

    :param graph: (networkx.Graph)
    :param node_index: (int)
    :return: (list(int))
    """
    for i, j, is_add_edge in graph.edges(data='add'):

        # No need to consider the reverse as both
        if is_add_edge and i == node_index:
            yield j

    return StopIteration


def add_reaction_edges(graph, reaction):
    """
    Add edges between reactants and products on a graph. For example a
    reversible reaction R <-> P would add two directional edges

    :param graph: (networkx.Graph)
    :param reaction: (rksim.systems.Reaction)
    """

    is_reversible = reaction.is_reversible()

    for reactant in reaction.reactants():
        for product in reaction.products():

            r_idx = name_to_index(reactant.name, graph)
            p_idx = name_to_index(product.name, graph)

            graph.add_edge(r_idx, p_idx, k=1.0, add=False)

            # If the reaction is reversible add a reaction edge from
            # each product to each reactant
            if is_reversible:
                graph.add_edge(p_idx, r_idx, k=1.0, add=False)

    return None


def add_addition_edges(graph, species_generator):
    """
    Add edges denoting an addition e.g. A + B -> P would have and addition edge
    between A and B

    :param graph: (networkx.Graph)
    :param species_generator: (generator)
    """

    species = list(species_generator)

    if len(species) == 2:
        # Add an addition edge forwards and backwards between the group
        for s1, s2 in (species, reversed(species)):
            graph.add_edge(name_to_index(s1.name, graph),
                           name_to_index(s2.name, graph),
                           add=True, k=None)

    if len(species) > 2:
        raise NotImplementedError

    return None


def name_to_index(name, graph):
    """Get the index of a node given it's name"""

    for i in graph.nodes:
        if graph.nodes[i]['name'] == name:
            return i

    return None


def add_nodes(graph, reactions):
    """
    Add nodes to a graph as species from a list of reactions

    :param graph: (networkx.Graph)
    :param reactions: (list(rksim.systems.Reaction))
    """
    # List to ensure no species is added more than once to the network
    added_names = []

    # Number of added nodes (iterator)
    n = 0

    for reaction in reactions:
        for species in reaction.components:

            if species.order != 1:
                # TODO allow larger orders
                raise NotImplementedError

            # Skip adding the node if already added
            if species.name in added_names:
                continue

            # Add this name to the list so it won't be added twice
            added_names.append(species.name)

            # Add the node and step the iterator
            graph.add_node(n,                       # index of the node
                           name=species.name,       # name of the node
                           c0=1E-8,                 # initial concentration
                           species=species)         # rksim.species.Species
            n += 1

    return None
