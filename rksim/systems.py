import numpy as np
from rksim.data import TimeSeries
import rksim.networks as nws
from rksim.plotting import plot


class System:

    def __str__(self):
        return f'{[species.name for species in self.species()]}'

    def set_rate_constants(self, k=1.0):
        """Set initial rate constants"""
        for edge in self.network.edges:
            self.network.edges[edge]['k'] = k

        return None

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
            species.simulated_time_series = TimeSeries(name=species.name,
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

        expt_series = [species.time_series for species in self.species()]
        sim_series = [species.simulated_time_series
                      for species in self.species()]

        return plot(expt_series, sim_series, name=name, dpi=dpi)

    def __init__(self, *args):
        """
        System of reactions

        :param args: (rksim.system.Reaction)
        """

        self.network = nws.make_network(args)
