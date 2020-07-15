from rksim.species import Reactant, Product


class Reaction:

    def _components(self, species_class):
        """Get the components of a particular class"""

        for species in self.components:
            if isinstance(species, species_class):
                yield species

        return None

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

        # Sort the components alphabetically by name
        self.components = sorted(self.components, key=lambda s: s.name)
        return None

    def __init__(self, *args):
        """Reaction e.g. R -> P"""

        self.components = []
        self.set_components(args)


class Irreversible(Reaction):
    """Irreversible reaction"""


class Reversible(Reaction):
    """Reversible reaction"""
