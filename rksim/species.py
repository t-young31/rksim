class Species:

    def __str__(self):
        return self.name

    def __eq__(self, other):
        return self.name == other.name

    def __init__(self, name=None):
        """
        Generic molecular species e.g. both R and P in R -> P

        :param name: (str) Name of this species
        """
        self.name = name

        # Stoichiometry of species in a reaction e.g. R + R -> P,
        # R.stoichiometry = 2 and P.stoichiometry = 1
        self.stoichiometry = 1

        self.series = None                  # rksim.data.TimeSeries
        self.simulated_series = None        # rksim.data.TimeSeries


class Reactant(Species):
    """Reactant species e.g. R in R -> P"""


class Product(Species):
    """Product species e.g. P in R -> P"""
