import numpy as np
from rksim.data import TimeSeries
from rksim.plotting import plot


class Species:

    def __init__(self, name=None):
        """
        Generic molecular species e.g. both R and P in R -> P

        :param name: (str) Name of this species
        """
        self.name = name

        self.time_series = None                 # rksim.data.TimeSeries
        self.simulated_time_series = None       # rksim.data.TimeSeries


class Reactant(Species):
    """Reactant species e.g. R in R -> P"""


class Product(Species):
    """Product species e.g. P in R -> P"""


class Reaction:

    def __init__(self, *args):
        """Reaction e.g. R -> P"""

        self.reactants = [mol for mol in args if isinstance(mol, Reactant)]
        self.products = [mol for mol in args if isinstance(mol, Product)]


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
            for species in reaction.reactants + reaction.products:
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
