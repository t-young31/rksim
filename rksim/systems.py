import numpy as np
from rksim.data import TimeSeries
import rksim.networks as nws
from rksim.plotting import plot
from rksim.exceptions import CannotSetAttribute
from rksim.exceptions import CannotGetAttribute


class System:

    def __str__(self):
        return f'{[species.name for species in self.species()]}'

    def mse(self):
        """
        Calculate the mean squared error between the
        true concentrations and the simulated concentrations

        :return: (float)
        """
        mse = 0

        for species in self.species():

            # Only compute the error on species with a time series
            if species.series is None:
                continue

            # For each time compute the difference between the simulated
            # concentration and the actual concentration. The simulated
            # concentration need to be the one closet to the current t
            # as the real and simulated series could be over different times
            for i, time in enumerate(species.series.times):

                idx = (np.abs(species.simulated_series.times - time)).argmin()

                diff = (species.series.concentrations[i] -
                        species.simulated_series.concentrations[idx])

                # Compute the square of the difference to make it smooth
                mse += diff**2

        return mse

    def rate_constants(self):
        """Get a numpy array of rate constants"""
        ks = []

        for edges in self.network.edge_mapping.values():
            # All the edges in this set have the same k, so use the first
            i, j = next(iter(edges))
            ks.append(self.network[i][j]['k'])

        return np.array(ks)

    def set_rate_constants(self, ks=None, k=None):
        """Set rate constants

        :param ks: (np.ndarray) shape = (n,) where n is the number of unique
                   rate constants in the reaction network

        :param k: (float) Rate constant to set all ks as
        """
        # Number of unique edges
        n_edges = len(self.network.edge_mapping)

        ks = ks if ks is not None else n_edges * [k]
        assert len(ks) == n_edges

        # Set all the equivalent edges with the specified rate constants
        for n, edges in enumerate(self.network.edge_mapping.values()):
            for (i, j) in edges:
                self.network[i][j]['k'] = ks[n]

        return None

    def set_initial_concentration(self, name, c):
        """
        Set the initial concentration (c0) for a species, given
        its name. Will set the c0 node attribute

        :param name: (str) Name of the species
        :param c: (float) Concentration (mol dm^-3)
        """
        try:
            node_index = self.network.node_mapping[name]

        except KeyError:
            raise CannotSetAttribute('Species not found in the network')

        # Set the value
        self.network.nodes[node_index]['c0'] = float(c)
        return None

    def rate_constant(self, name_i, name_j):
        """
        Set the rate constant (k) between a species i and j in the
        reaction network is directional so i -> j rate constant

        :param name_i: (str) Name of a species
        :param name_j: (str) Name of a species
        """
        try:
            edge = self.network.get_edge(name_i, name_j)
            return edge['k']

        except KeyError:
            raise CannotGetAttribute('Reaction not found in the network')

    def set_rate_constant(self, name_i, name_j, k):
        """
        Set the rate constant (k) between a species i and j in the
        reaction network is directional so i -> j rate constant

        :param name_i: (str) Name of a species
        :param name_j: (str) Name of a species
        :param k: (float) Rate constant
        """
        try:
            edge = self.network.get_edge(name_i, name_j)
            edge['k'] = k

        except KeyError:
            raise CannotSetAttribute('Reaction not found in the network')

    def set_simulated(self, concentrations, times):
        """
        Set concentrations as a function of time for all species in this
        system of reactions

        :param concentrations: (np.ndarray) Array of concentrations (mol dm^-3)
                               for all components in this system. shape (n, m)
                               where there are m components in this reaction

        :param times: (np.ndarray) Array of times in s. shape = (n,)
        """
        for i, species in enumerate(self.species()):

            # Concentrations are a matrix of time points as the rows and
            # columns as the different species
            concs = concentrations[:, i]
            species.simulated_series = TimeSeries(name=species.name,
                                                  times=times,
                                                  concentrations=concs)
        return None

    def derivative(self, concentrations, time=0.0):
        """
        Calculate the derivative of all the concentrations with respect to time

        :param time: (float) Time in s at which the derivative is calculated
                     ((needs to be a parameter for odeint))

        :param concentrations: (np.ndarray) array of concentrations in mol
                               dm^-3 shape = (n,) where n is the number of
                               components (species in this system). Must be >0
        """
        n = len(concentrations)
        dcdt = np.zeros(n)

        for i in range(n):
            inflows = 0
            outflows = 0

            neighbours_i = list(self.network.neighbours_a(i))

            # Add the rate constant for all outflowing reactions
            #          P1
            #          ^
            #          |
            # i.e. A + B -> P2    then outflows for B is (k1 + k2)
            # and j in [P1, P2]
            for j in nws.outflow_node(self.network, i):
                # Add the rate constant for this outflow reaction
                # will be multiplied by [A][B] and divided by the number
                # of neighbours + 1 i.e. the above n_neighbours = 0

                # Stoichiometry from i -> j
                sto = self.network.edges[(i, j)]['sto'][0]
                conc = sto * concentrations[i]**sto

                # Concentrations of components on addition edges to this
                # node that has a reaction edge to j need to be multiplied
                for k in neighbours_i:

                    if (k, j) not in self.network.edges:
                        continue

                    # Only add reaction edges between neighbours to j and i
                    if self.network.edges[(k, j)]['add'] is True:
                        continue

                    sto = self.network.edges[(k, j)]['sto'][0]
                    conc *= sto * concentrations[k]**sto

                # Number of nodes that this outflow is distributed over
                m = 1
                for k in self.network.neighbours_a(j):
                    if (i, k) in self.network.edges:
                        m += 1

                # Add the outflow as k[A]^o[B]^p
                outflows += self.network[i][j]['k'] * conc / m

            # Add the rate constant for all inflowing reactions to this
            # node for example:  A + B -> P then inflow nodes are A and B
            for j in nws.inflow_node(self.network, i):

                # Stoichiometry from j -> i as the second element in the tuple
                sto_j, sto_i = self.network.edges[(j, i)]['sto']
                conc = sto_i * concentrations[j]**sto_j

                # Number of nodes that this inflow comes from
                m = 1
                for k in self.network.neighbours_a(j):

                    if (k, i) not in self.network.edges:
                        continue

                    # Only add reaction edges between neighbours to k and i
                    if self.network.edges[(k, i)]['add'] is True:
                        continue

                    m += 1
                    sto_k, sto_i = self.network.edges[(k, i)]['sto']
                    conc *= sto_i * concentrations[k]**sto_k

                inflows += conc * self.network[j][i]['k'] / m

            # Derivative is the inflow minus the outflow e.g.
            # R -> P  d[R]/dt = -k[R], d[P]/dt = k[R]
            dcdt[i] = inflows - outflows

        return dcdt

    def species(self):
        """Get the next species in this system from the reaction network"""

        for i in self.network.nodes:
            yield self.network.nodes[i]['species']

        return None

    def plot(self, name='system', dpi=400):
        """Plot both the simulated and experimental data for this system"""
        return plot(self.species(), name=name, dpi=dpi)

    def __init__(self, *args):
        """
        System of reactions

        :param args: (rksim.system.Reaction)
        """

        self.network = nws.make_network(args)
