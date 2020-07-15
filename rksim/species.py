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

        self.series = None                  # rksim.data.TimeSeries
        self.simulated_series = None        # rksim.data.TimeSeries


class Reactant(Species):
    """Reactant species e.g. R in R -> P"""


class Product(Species):
    """Product species e.g. P in R -> P"""
