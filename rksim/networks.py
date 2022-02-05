import networkx as nx
from rksim.plotting import show_or_plot
from copy import deepcopy


class Network(nx.DiGraph):

    def add_nodes(self, reactions):
        """
        Add nodes to a graph as species from a list of reactions

        :param self: (rksim.networks.Network)
        :param reactions: (list(rksim.systems.Reaction))
        """
        # List to ensure no species is added more than once to the network
        added_names = []

        # Number of added nodes (iterator)
        n = 0

        for reaction in reactions:
            for species in deepcopy(reaction.components):

                # A species only has a stoichiometry in a reaction
                species.stoichiometry = None

                # Skip adding the node if already added
                if species.name in added_names:
                    continue

                # Add this name to the list so it won't be added twice
                added_names.append(species.name)

                # Add the node and step the iterator
                self.add_node(n,  # index of the node
                              name=species.name,  # name of the node
                              c0=1E-15,           # initial concentration
                              species=species)    # rksim.species.Species
                n += 1

        # Set the mapping between node names and node indexes
        self.set_node_mapping()
        return None

    def set_node_mapping(self):
        """Set the mapping from node names to indexes"""
        for i in self.nodes:
            name = self.nodes[i]['name']
            self.node_mapping[name] = i

        return None

    def plot(self, name=None):
        """Plot a reaction network using NetworkX and matplotlib"""
        nx.draw_networkx(self)

        return show_or_plot(name)

    def __init__(self, *args):
        """Subclass of networkx.Graph. Once initialised then calling
         self.add_node() will break the mappings"""
        super().__init__()

        self.node_mapping = {}              # Mapping from names -> indexes
        self.add_nodes(args)
