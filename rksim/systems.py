import numpy as np
from rksim.data import TimeSeries
import rksim.networks as nws
from rksim.reactions import ReactionSet
from rksim.plotting import plot
from rksim.exceptions import CannotSetAttribute
from rksim.exceptions import CannotGetAttribute


class System:

    def __str__(self):
        return f'{[species.name for species in self.species()]}'

    def mse(self, relative=False):
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

                # Compute the square of the difference to make the function
                # more smooth. Compute the relative error for a more
                # reasonable fit to low and high concentration data
                if relative and species.series.concentrations[i] > 0:
                    mse += diff ** 2 / species.series.concentrations[i]

                else:
                    mse += diff ** 2

        return mse

    def rate_constants(self):
        """Get a numpy array of rate constants"""
        return np.array([reaction.k for reaction in self.reactions])

    def rate_constant(self, *args):
        """Get a rate constant for a reaction.

        system.rate_constant('R', 'P') :-> k_RP if the reaction R -> P is in
        the system

        :param args: (str) >1 Name of a species in this system
        """
        return self.reactions.rate_constant(species_names=args)

    def set_rate_constant(self, *args, k=1.0):
        """Set a rate constant for a reaction"""
        return self.reactions.set_rate_constant(species_names=args, k=k)

    def set_rate_constants(self, ks=None, k=None):
        """Set rate constants

        :param ks: (np.ndarray) shape = (n,) where n is the number of unique
                   rate constants in the reaction network

        :param k: (float) Rate constant to set all ks as
        """
        return self.reactions.set_rate_constants(ks=ks, k=k)

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
        for i, species in enumerate(self.species()):

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
                                                      concentrations=concs)
        return None

    def set_stos(self):
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

    def derivative(self, concentrations, time=0.0):
        """
        Calculate the derivative of all the concentrations with respect to time

        :param time: (float) Time in s at which the derivative is calculated
                     ((needs to be a parameter for odeint))

        :param concentrations: (np.ndarray) array of concentrations in mol
                               dm^-3 shape = (n,) where n is the number of
                               components (species in this system). Must be >0
        """
        # TODO optimise this function
        n = len(concentrations)
        dcdt = np.zeros(n)

        for i in range(n):
            inflow, outflow = 0.0, 0.0

            for j, reaction in enumerate(self.reactions):

                # If this component is outflowing
                if self.stos[j, i, 0] != 0:
                    sto = self.stos[j, i, 0]
                    conc = sto * concentrations[i] ** sto

                    for k in range(n):

                        # Don't re-multiply by this conc and only multiply if
                        # the stoichometry of a reactant is > 0
                        if k == i or self.stos[j, k, 0] == 0:
                            continue

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

                        # Reactants flowing to i
                        if self.stos[j, k, 0] == 0:
                            continue

                        conc *= concentrations[k] ** self.stos[j, k, 0]

                    inflow += reaction.k * conc

            # Set the derivative dC_i/dt
            dcdt[i] = inflow - outflow

        return dcdt

    def species(self):
        """Get the next species in this system from the reaction network"""

        for i in self.network.nodes:
            yield self.network.nodes[i]['species']

        return None

    def fit(self, data, optimise=True):
        """Fit some data to this system. i.e. optimise ks"""
        return data.fit(self, optimise=optimise)

    def plot(self, name='system', dpi=400):
        """Plot both the simulated and experimental data for this system"""
        return plot(self.species(), name=name, dpi=dpi)

    def __init__(self, *args):
        """
        System of reactions

        :param args: (rksim.system.Reaction)
        """

        self.network = nws.Network(*args)
        self.reactions = ReactionSet(*args)

        self.stos = None
        self.set_stos()
