import networkx as nx
from rksim.plotting import show_or_plot


class Network(nx.DiGraph):

    def add_addition_edges(self, species):
        """
        Add edges denoting an addition e.g. A + B -> P would have and addition
        edge between A and B

        :param species: (generator)
        """
        species = list(species)

        if len(species) == 2:
            # Add an addition edge forwards and backwards between the group
            for s1, s2 in (species, reversed(species)):
                self.add_edge(node_name_to_index(s1.name, self),
                              node_name_to_index(s2.name, self),
                              add=True, k=None)

        if len(species) > 2:
            raise NotImplementedError

        return None

    def add_reaction_edges(self, reaction, default_k=1.0):
        """
        Add edges between reactants and products on a graph. For example a
        reversible reaction R <-> P would add two directional edges

        :param reaction: (rksim.systems.Reaction)
        :param default_k: (float) Default rate constant for this reaction
        """
        is_reversible = reaction.is_reversible()

        for reactant in reaction.reactants():
            for product in reaction.products():

                r_idx = node_name_to_index(reactant.name, self)
                p_idx = node_name_to_index(product.name, self)

                self.add_edge(r_idx, p_idx, k=default_k, add=False)

                # If the reaction is reversible add a reaction edge from
                # each product to each reactant
                if is_reversible:
                    self.add_edge(p_idx, r_idx, k=default_k, add=False)

        return None

    def set_edge_mapping(self):
        """
        Set the mapping between unique edges and all edges. For example,
        for a A + B -> P reaction with a reaction network

              0
            A --\
        1  +|     P
            B - /
              2

        then the edges A-P and B-P need to have identical rate constants (k).
        This function will generate a dictionary keyed with one edge (e.g. A-P)
        and valued with all the equivalent edges. For the above this would be
        {0: [(A, P), (B, P)]} where 0 is the index of the A-P edge and (A, P)
        are actually the node numbers in the graph
        """
        self.edge_mapping = {}
        added_edges = []

        for n, edge in enumerate(self.edges):
            i, j = edge

            # Only consider reaction edges
            if self.edges[edge]['add'] is True:
                continue

            # Don't add this edge if it's already in the mapping
            if edge in added_edges:
                continue

            # At least this edge has the same rate constant
            edges_same_k = [edge]

            # Iterate through the neighbours to node i, above if i = A then
            # B is a neighbour
            for neighbour in neighbours(self, i):

                # If the neighbour has a reaction edge to the same end point
                # (which cannot be an add edge) then the rate constants need
                # to be the same for this current edge and (neighbour, j)
                if (neighbour, j) in self.edges:
                    edges_same_k.append((neighbour, j))

            # Add reaction edges that are neighbours to both i and j e.g.
            # A + B -> C + D where i = A and j = C then this adds an edge
            # to edges_same_k between B and D
            for neighbour_i in neighbours(self, i):
                for neighbour_j in neighbours(self, j):
                    n_edge = (neighbour_i, neighbour_j)

                    # If this reaction edge exists
                    if n_edge in self.edges and not self.edges[n_edge]['add']:
                        edges_same_k.append(n_edge)

            # Likewise with neighbours to j if i = R for a reaction R -> C + D
            # here i = R and j = C then an edges_same_k will be added between
            # R and D
            for neighbour in neighbours(self, j):
                if (i, neighbour) in self.edges:
                    edges_same_k.append((i, neighbour))

            # Set the mapping dictionary
            self.edge_mapping[n] = edges_same_k
            # and update which edges have been added so far
            added_edges += edges_same_k

        return None

    def plot(self, name=None, dpi=400):
        """Plot a reaction network using NetworkX and matplotlib"""
        nx.draw_networkx(self)

        return show_or_plot(name, dpi)

    def __init__(self):
        """Subclass of networkx.Graph"""
        super().__init__()

        self.edge_mapping = None


def make_network(reactions):
    """
    Make a reaction network for a system of reactions. Will have 'add' edges
    and 'reaction' edges, where the former are denoted with '+'

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
    network = Network()

    # Add all species as nodes to the graph/network
    add_nodes(network, reactions=reactions)

    for reaction in reactions:
        # Add any + edges between reactants, then products
        network.add_addition_edges(reaction.reactants())
        network.add_addition_edges(reaction.products())

        # Add reaction arrow(s) with associated rate constant (k) values
        network.add_reaction_edges(reaction)

    # Set the number of + edges on each node
    populate_neighbour_number(network)

    # Set the mapping for unique rate constants onto edges
    network.set_edge_mapping()

    return network


def inflow_node(network, node_index):
    """
    Get a nodes that flows into a node_index

    :param network: (networkx.Graph)
    :param node_index: (int)
    """
    for j in network.predecessors(node_index):
        if network[j][node_index]['add'] is False:
            yield j

    return StopIteration


def outflow_node(network, node_index):
    """
    Get a nodes that flows from a node_index

    :param network: (networkx.Graph)
    :param node_index: (int)
    """
    for j in network.successors(node_index):
        if network[node_index][j]['add'] is False:
            yield j

    return StopIteration


def populate_neighbour_number(network):
    """
    Populate the number of neighbours over add edges for a reaction graph

           P
          ^ ^
         /  \
       /     \
      A ----- B
          +

    P_n_neighbours = 0, A_n_neighbours = 1, B_n_neighbours = 1

    :param network: (networkx.Graph)
    """
    neighbours_dict = {}

    for node in network.nodes:
        n_neighbours = 0

        for neighbour in network.neighbors(node):
            if network[node][neighbour]['add']:
                n_neighbours += 1

        neighbours_dict[node] = n_neighbours

    # Set the number of neighbours for each species in the network/graph
    nx.set_node_attributes(network, neighbours_dict, name='n_neighbours')
    return None


def neighbours(network, node_index):
    """
    Get the neighbours along addition edges next to this atom. For example
    a graph

        +
    A ----- B

    where A is indexed 0 and B 1 then get_add_neighbours(graph, 0) -> [1]

    :param network: (networkx.Graph)
    :param node_index: (int)
    :return: (list(int))
    """
    for i, j, is_add_edge in network.edges(data='add'):

        # No need to consider the reverse as both
        if is_add_edge and i == node_index:
            yield j

    return StopIteration


def node_name_to_index(name, graph):
    """Get the index of a node given it's name"""

    for i in graph.nodes:
        if graph.nodes[i]['name'] == name:
            return i

    return None


def add_nodes(network, reactions):
    """
    Add nodes to a graph as species from a list of reactions

    :param network: (networkx.Graph)
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
            network.add_node(n,  # index of the node
                             name=species.name,  # name of the node
                             c0=1E-8,  # initial concentration
                             species=species)         # rksim.species.Species
            n += 1

    return None
