from rksim.species import Reactant, Product
import rksim.exceptions as ex
from copy import deepcopy


def prune_species(components):
    """Remove any duplicates and set stoichiometries"""
    unique_names = set(species.name for species in components)
    unique_species = []

    for name in unique_names:
        identical_species = [s for s in components if s.name == name]

        # Add this species only once
        species = identical_species[0]

        # The order in this species is the number of times it appears
        # as either a reactant or product
        species.stoichiometry = len(identical_species)

        unique_species.append(species)

    # Sort the components alphabetically by name
    return sorted(unique_species, key=lambda s: s.name)


class Reaction:

    def _components(self, species_class):
        """Get the components of a particular class"""

        for species in self.components:
            if isinstance(species, species_class):
                yield species

        return StopIteration

    def swap_reactants_and_products(self):
        """Switch reactants to products e.g. R -> P :-> P -> R"""
        for species in self.components:
            if isinstance(species, Reactant):
                species.__class__ = Product

            elif isinstance(species, Product):
                species.__class__ = Reactant

            else:
                raise ex.RKSimCritical('Had unknown species in reaction')

        return None

    def ordered_names(self):
        reactant_names = [r.name for r in self.reactants()]
        product_names = [p.name for p in self.products()]
        return reactant_names + product_names

    def is_reversible(self):
        """Is this reaction reversible?"""

        if isinstance(self, Reversible):
            return True

        if isinstance(self, Irreversible):
            return False

    def reactants(self):
        """Get the next reactant in this reaction"""
        return self._components(Reactant)

    def products(self):
        """Get the next reactant in this reaction"""
        return self._components(Product)

    def sto(self, name):

        for species in self.components:
            if species.name == name:
                return species.stoichiometry

        # If this
        return 0

    def __init__(self, *args):
        """Reaction e.g. R -> P"""
        self.components = prune_species(args)
        self.k = 1.0                            # Rate constant


class Irreversible(Reaction):
    """Irreversible reaction"""


class Reversible(Reaction):
    """Reversible reaction"""


class ReactionSet:

    def __iter__(self):
        return iter(self._list)

    def __len__(self):
        return len(self._list)

    def species(self):

        for reaction in self._list:
            for species in reaction.components:
                yield species

        return StopIteration

    def rate_constant(self, species_names):
        """Get a named rate constant"""

        for reaction in self._list:
            if list(species_names) == list(reaction.ordered_names()):
                return reaction.k

        raise ex.CannotGetAttribute('Reaction not found')

    def set_rate_constant(self, species_names, k=1.0):
        """Get a named rate constant"""

        for reaction in self._list:
            if list(species_names) == list(reaction.ordered_names()):
                reaction.k = k
                return

        raise ex.CannotSetAttribute('Reaction not found')

    def set_rate_constants(self, ks=None, k=None):
        """Set the rate constants either with an ordered list of ks, or
        all the same rate constant

        :param ks: (list(float)) or (np.ndarray) shape (len(self.reactions),)
        :param k: (float)
        """
        n = len(self._list)

        ks = n * [k] if k is not None else ks
        if len(ks) != n:
            raise ex.CannotSetAttribute('Incorrect number of rate constants')

        # Set the rate constants for all reactions
        for i, reaction in enumerate(self._list):
            reaction.k = ks[i]

        return None

    def rate_constants(self):
        """Get all the rate constants"""

        for reaction in self._list:
            yield reaction.k

        return StopIteration

    def add_reverse_reactions(self):
        """If there are reversible reactions swap for irreversible
        and add the corresponding reverse reaction"""
        reactions = []

        for reaction in self._list:

            # Don't modify irreversible forward reactions
            if not reaction.is_reversible():
                reactions.append(reaction)
                continue

            # Add the forward reaction
            reactions.append(reaction)

            # Copy the forward reaction, swap it and append to the list
            swapped_reaction = deepcopy(reaction)
            swapped_reaction.swap_reactants_and_products()

            reactions.append(swapped_reaction)

        # Reset the reactions with the fully populated list
        self._list = reactions
        return None

    def __init__(self, *args):
        """Set of reactions"""
        self._list = args
        self.add_reverse_reactions()
