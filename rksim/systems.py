import numpy as np
from rksim.data import TimeSeries
import rksim.networks as nws
from rksim.fit import fit
from rksim.reactions import ReactionSet
from rksim.plotting import plot
from rksim.exceptions import CannotSetAttribute


class System(ReactionSet):

    def rmse(self, relative=False):
        """Calculate the root mean square error over the tim series"""
        return np.sqrt(self.mse(relative=relative))

    def mse(self, relative=False):
        """
        Calculate the mean squared error between the
        true concentrations and the simulated concentrations

        :return: (float)
        """
        mse = 0

        for species in self.species:

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

                # Compute the square of the difference to make the function
                # more smooth. Compute the relative error for a more
                # reasonable fit to low and high concentration data
                if relative and species.series.concentrations[i] > 0:
                    mse += diff ** 2 / species.series.concentrations[i]

                else:
                    mse += diff ** 2

        return mse

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

    def set_simulated(self, concentrations, times):
        """
        Set concentrations as a function of time for all species in this
        system of reactions

        :param concentrations: (np.ndarray) Array of concentrations (mol dm^-3)
                               for all components in this system. shape (n, m)
                               where there are m components in this reaction

        :param times: (np.ndarray) Array of times in s. shape = (n,)
        """
        for i, species in enumerate(self.species):

            # Concentrations are a matrix of time points as the rows and
            # columns as the different species
            concs = concentrations[:, i]

            # If a time series is already set only update the concentrations
            if (species.simulated_series is not None
                    and len(concs) == len(species.simulated_series.times)):
                species.simulated_series.concentrations = concs

            # Otherwise set the time series
            else:
                species.simulated_series = TimeSeries(name=species.name,
                                                      times=times,
                                                      concentrations=concs,
                                                      simulated=True)
        return None

    def set_stoichiometries(self):
        """Set an array of stoichiometries"""
        n = len(self.reactions)
        m = self.network.number_of_nodes()

        # Matrix of stoichiometries (sto). Reactions as rows and unique
        # species as columns
        self.stos = np.zeros(shape=(n, m, 2))

        for i, reaction in enumerate(self.reactions):
            reactant_names = [r.name for r in reaction.reactants()]

            for species_name, j in self.network.node_mapping.items():

                # Final idx is 0 if this species is a Reactant or 1 if Product
                k = 0 if species_name in reactant_names else 1

                self.stos[i, j, k] = reaction.sto(species_name)

        return None

    def component_derivative(self, i, concentrations):
        """Calculate the derivative with respect to a component i in the
        system.
                                      __
        dc_i/dt = (   Σ   S_{j,i} k_j  || c_k ^S_{j, k}  -
                  ( j ∈ P              k
                                      __
                     Σ   S_{j,i} k_j  || c_k ^S_{j, k} )
                    j ∈ R              k               )

        where P is the set of reactions where i is a product and R is the
        set of reactions where i is a reactant
        """
        n = len(concentrations)
        inflow, outflow = 0.0, 0.0

        for j, reaction in enumerate(self.reactions):

            # If this component is outflowing
            if self.stos[j, i, 0] != 0:
                conc = self.stos[j, i, 0]

                for k in range(n):

                    # If the reactant is not present then the stoichometry
                    # is 0 and conc will be multipled by 1
                    conc *= concentrations[k] ** self.stos[j, k, 0]

                # Add e.g. k[A][B]^2 to the outflow for this reaction
                outflow += reaction.k * conc

            # If this component is inflow
            if self.stos[j, i, 1] != 0:

                # Product concentration is just the stoichiometry e.g. for
                # i == C in A -> 2C then sto = 2
                conc = self.stos[j, i, 1]

                # For all the reactants in this reaction multiply by
                # their concentration raised to the power of their
                # stoichiometry
                for k in range(n):
                    conc *= concentrations[k] ** self.stos[j, k, 0]

                inflow += reaction.k * conc

        return inflow - outflow

    def derivative(self, concentrations, time=0.0):
        """
        Calculate the derivative of all the concentrations with respect to time

        :param time: (float) Time in s at which the derivative is calculated
                     ((needs to be a parameter for odeint))

        :param concentrations: (np.ndarray) array of concentrations in mol
                               dm^-3 shape = (n,) where n is the number of
                               components (species in this system). Must be >0
        """
        n_components = len(concentrations)
        return np.array([self.component_derivative(i, concentrations)
                         for i in range(n_components)])

    @property
    def species(self):
        """Get the next species in this system from the reaction network"""

        for i in self.network.nodes:
            yield self.network.nodes[i]['species']

        return None

    def simulate(self, max_time):
        """Simulate this system to time = max_time"""

        if all(self.network.nodes[i]['c0'] < 1E-8 for i in self.network.nodes):
            raise RuntimeError('Cannot simulate a system with all zero '
                               'concentrations. Set some concentrations with'
                               ' set_initial_concentration()')

        return self.fit(data=None, optimise=False, max_time=max_time)

    def fit(self, data, optimise=True, max_time=None):
        """Fit some data to this system. i.e. optimise ks"""
        if data is not None:
            data.assign(self)

        return fit(data, self, optimise, max_time)

    def plot(self, name='system',exc_species=None):
        """Plot both the simulated and experimental data for this system"""
        species = self.species

        if exc_species is not None:
            species = [s for s in species if s.name not in exc_species]

        return plot(species, name=name)

    def __init__(self, *args):
        """
        System of reactions. Subclass of ReactionSet with a self.reactions
        attribute

        :param args: (rksim.system.Reaction)
        """
        super().__init__(*args)

        # Network of unique components in the system
        self.network = nws.Network(*args)

        # Stoichiometry matrix (np.ndarray)
        self.stos = None
        self.set_stoichiometries()
