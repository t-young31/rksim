import numpy as np
from rksim.data import TimeSeries
import rksim.networks as nws
from rksim.plotting import plot
from rksim.exceptions import CannotSetAttribute


class System:

    def __str__(self):
        return f'{[species.name for species in self.species()]}'

    def get_mse(self):
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

    def get_rate_constants(self):
        """Get a numpy array of rate constants"""

        # Make a mapping....
        raise NotImplementedError

        return

    def set_rate_constants(self, k=1.0):
        """Set initial rate constants"""

        # self.network.edges[edge]['k'] = k

        return None

    def set_initial_concentration(self, name, c):
        """
        Set the initial concentration (c0) for a species, given
        its name. Will set the c0 node attribute

        :param name: (str) Name of the species
        :param c: (float) Concentration (mol dm^-3)
        """

        for i in self.network.nodes:
            node = self.network.nodes[i]

            if node['name'] == name:
                node['c0'] = float(c)
                return

        raise CannotSetAttribute('Species not found in the network')

    def set_rate_constant(self, name_i, name_j, k):
        """
        Set the rate constant (k) between a species i and j in the
        reaction network is directional so i -> j rate constant

        :param name_i: (str) Name of a species
        :param name_j: (str) Name of a species
        :param k: (float) Rate constant
        """

        for edge in self.network.edges:

            i, j = edge
            # If the names of the nodes are those specified then
            # set the rate constant between them
            if (self.network.nodes[i]['name'] == name_i and
                    self.network.nodes[j]['name'] == name_j):

                self.network.edges[edge]['k'] = k
                return

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

            concs = concentrations[:, i]
            species.simulated_series = TimeSeries(name=species.name,
                                                  times=times,
                                                  concentrations=concs)
        return None

    def derivative(self, concentrations, time=0.0):
        """
        Calculate the derivative of all the concentrations with respect to time

        :param time: (float) Time in s at which the derivative is calculated

        :param concentrations: (np.ndarray) array of concentrations in mol
                               dm^-3 shape = (n,) where n is the number of
                               components (species in this system). Must be >0
        """
        n = len(concentrations)
        dcdt = np.zeros(n)

        for i in range(n):
            inflows = 0
            outflows = 0

            # Concentration on this node
            conc = concentrations[i]

            # Product of concentrations on this node e.g. [A][B]
            # if i == A and the reaction is A + B -> P
            for neighbour in nws.neighbours(self.network, i):
                conc *= concentrations[neighbour]

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
                n_neighbours = self.network.nodes[j]['n_neighbours']

                outflows += self.network[i][j]['k'] / (n_neighbours + 1)

            outflows *= conc

            # Add the rate constant for all inflowing reactions to this
            # node for example:  A + B -> P then inflow is
            for j in nws.inflow_node(self.network, i):

                n_j = self.network.nodes[j]['n_neighbours']
                rate_constant = self.network[j][i]['k']

                inflow = concentrations[j] * rate_constant / (n_j + 1)

                for k in nws.neighbours(self.network, j):
                    inflow *= concentrations[k]

                inflows += inflow

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

        expt_series = [species.series for species in self.species()]
        sim_series = [species.simulated_series for species in self.species()]

        return plot(expt_series, sim_series, name=name, dpi=dpi)

    def __init__(self, *args):
        """
        System of reactions

        :param args: (rksim.system.Reaction)
        """

        self.network = nws.make_network(args)
