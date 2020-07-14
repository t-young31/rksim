import numpy as np
from rksim.data import TimeSeries
from rksim.graphs import make_network
from rksim.plotting import plot


class Species:

    def __eq__(self, other):
        return self.name == other.name

    def __init__(self, name=None):
        """
        Generic molecular species e.g. both R and P in R -> P

        :param name: (str) Name of this species
        """
        self.name = name

        # Order of species in a reaction e.g. R + R -> P, R.order = 2
        # and P.order = 1
        self.order = 1

        self.time_series = None                 # rksim.data.TimeSeries
        self.simulated_time_series = None       # rksim.data.TimeSeries


class Reactant(Species):
    """Reactant species e.g. R in R -> P"""


class Product(Species):
    """Product species e.g. P in R -> P"""


class Reaction:

    def set_components(self, components):
        """Remove any duplicates and set orders"""
        unique_names = set(species.name for species in components)

        for name in unique_names:
            identical_species = [s for s in components if s.name == name]

            # Add this species only once
            species = identical_species[0]

            # The order in this species is the number of times it appears
            # as either a reactant or product
            species.order = len(identical_species)

            self.components.append(species)

        return None

    def __init__(self, *args):
        """Reaction e.g. R -> P"""

        self.components = []
        self.set_components(args)


class Irreversible(Reaction):
    kf = 1.0


class Reversible(Reaction):
    kf = 1.0
    kb = 1.0


class System:

    def set_init_ks(self):
        """Set initial rate constants"""
        # TODO something sensible here
        return None

    def set_simulated(self, concentrations, times):
        """
        Set concentrations as a function of time for

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

    def derivative(self, concentrations, t):
        """Calculate the derivative of all the concentrations wrt time"""
        if len(self.reactions) > 1:
            raise NotImplementedError

        rconc, p_conc = concentrations
        reaction = self.reactions[0]

        # TODO work out the general formulation of this
        assert isinstance(reaction, Irreversible)

        deriv = [-reaction.kf * rconc, reaction.kf * rconc]

        return np.array(deriv)

    def species(self):
        """Get the next species in this system"""

        for reaction in self.reactions:
            for species in reaction.components:
                yield species

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
        self.reactions = args

        self.network = make_network(self)
