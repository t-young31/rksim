from rksim.species import Reactant, Product
import rksim.exceptions as ex
import numpy as np
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

    def __str__(self):
        return (f'Reaction({"+".join(r.name for r in self.reactants())} -> '
                f'{"+".join(p.name for p in self.products())})')

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
        """Get a list """
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
        """Get the stoichiometry of a named component in this reaction"""

        for species in self.components:
            if species.name == name:
                return species.stoichiometry

        return 0

    def _init_from_string(self, string):
        """Initialise a reaction from a string e.g. A+B->C"""
        l_components, r_components = string.split('->')
        components = [Reactant(name) for name in l_components.split('+')]
        components += [Product(name) for name in r_components.split('+')]

        self.components = prune_species(components)
        return None

    def __init__(self, *args):
        """Reaction e.g. R -> P"""
        if (len(args) == 1 and not
            (isinstance(args[0], Reactant) or isinstance(args[0], Product))):
            self._init_from_string(args[0])

        else:
            self.components = prune_species(args)

        self.k = 1.0                            # Rate constant


class Irreversible(Reaction):
    """Irreversible reaction"""


class Reversible(Reaction):
    """Reversible reaction"""


class ReactionSet:

    def __str__(self):
        rxn_str = "\n\t".join(str(rxn) for rxn in self.reactions)
        return f'Reactions({rxn_str})'

    def __iter__(self):
        return iter(self.reactions)

    def __len__(self):
        return len(self.reactions)

    @property
    def species(self):
        """Yield a species in this reaction set. May have repeats"""

        for reaction in self.reactions:
            for species in reaction.components:
                yield species

        return StopIteration

    def rate_constant(self, *args):
        """
        Get a named rate constant from the constituent species e.g. for a
        reaction set containing A + B -> C get the rate constant with
        rset.rate_constant(['A', 'B', 'C']). ['B', 'A', 'C'] *wont* work

        :param species_names: (list(str))
        """

        for reaction in self.reactions:
            if list(args) == reaction.ordered_names():
                return reaction.k

        raise ex.CannotGetAttribute('Reaction not found')

    def set_rate_constant(self, *args, k):
        """Get a named rate constant

        :param args: (str | int) either the name of the reaction or the
                     index in the system
        :param k: (float) Rate constant
        """

        for i, reaction in enumerate(self.reactions):
            if list(args) == reaction.ordered_names():
                reaction.k = k
                return

            try:
                if int(args[0]) == i:
                    reaction.k = k
                    return

            except ValueError:
                # argument not an integer
                continue

        raise ex.CannotSetAttribute(f'Reaction not found: {args}')

    def set_rate_constants(self, ks=None, k=None):
        """Set the rate constants either with an ordered list of ks, or
        all the same rate constant

        :param ks: (np.ndarray) shape = (n,) where n is the number of unique
                   rate constants in the reaction network
        :param k: (float)
        """
        n = len(self.reactions)

        ks = n * [k] if k is not None else ks
        if len(ks) != n:
            raise ex.CannotSetAttribute('Incorrect number of rate constants')

        # Set the rate constants for all reactions
        for i, reaction in enumerate(self.reactions):
            reaction.k = ks[i]

        return None

    def rate_constants(self):
        """YGet a numpy array of rate constants"""
        return np.array([reaction.k for reaction in self.reactions])

    def add_reverse_reactions(self):
        """If there are reversible reactions swap for irreversible
        and add the corresponding reverse reaction"""
        reactions = []

        for reaction in self.reactions:

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
        self.reactions = reactions
        return None

    def __init__(self, *args):
        """Set of reactions"""

        self.reactions = args
        self.add_reverse_reactions()
