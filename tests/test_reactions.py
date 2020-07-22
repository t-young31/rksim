from rksim.reactions import Reaction, ReactionSet
from rksim.reactions import Irreversible, Reversible
from rksim.species import Reactant, Product, Species
import rksim.exceptions as ex
import pytest


def test_reaction_init():

    r, p = Reactant('R'), Product('P')
    reaction = Reaction(r, p)

    assert r in reaction.components
    assert p in reaction.components

    assert len(reaction.components) == 2
    assert len(list(reaction.reactants())) == 1
    assert len(list(reaction.products())) == 1

    # Stoichiometries of R and P in this reaction
    assert reaction.sto('R') == 1
    assert reaction.sto('P') == 1

    # Names of the reaction components sorted by reactants -> products
    assert reaction.ordered_names() == ['R', 'P']

    irrev = Irreversible(r, p)
    assert not irrev.is_reversible()

    rev = Reversible(r, p)
    assert rev.is_reversible()


def test_swapping():

    r, p = Reactant('R'), Product('P')
    reaction = Reaction(r, p)
    reaction.swap_reactants_and_products()

    assert reaction.ordered_names() == ['P', 'R']

    # Cannot swap reactants and products if there are unknown species
    reaction2 = Reaction(r, Species('X'))
    with pytest.raises(ex.RKSimCritical):
        reaction2.swap_reactants_and_products()


def test_pruning():

    reaction = Irreversible(Reactant('A'), Reactant('A'), Product('P'))

    # Should remove the duplicate reactant
    assert len(reaction.components) == 2
    assert reaction.sto('A') == 2

    # Stoichiometry of an unknown component is 0
    assert reaction.sto('X') == 0


def test_reaction_set():

    rxn_set = ReactionSet(Reaction(Reactant('R'), Product('P')))

    assert len(rxn_set) == 1
    assert len(list(rxn_set.species())) == 2
    assert len(list(rxn_set.rate_constants())) == 1.0
    assert all(k == 1.0 for k in rxn_set.rate_constants())

    # Rate constants are initialised at 1
    assert rxn_set.rate_constant('R', 'P') == 1.0

    rxn_set2 = ReactionSet(Reversible(Reactant('R'), Product('P')))
    # Reaction set should add the reverse reaction
    assert len(rxn_set2) == 2

    rxn_set2.set_rate_constant('R', 'P', k=2.0)
    assert rxn_set2.rate_constant('R', 'P') == 2.0
    # Reverse reaction should still have k = 1.0
    assert rxn_set2.rate_constant('P', 'R') == 1.0

    # Cannot set a reaction set with a single reaction with two rate constants
    with pytest.raises(ex.CannotSetAttribute):
        rxn_set.set_rate_constants(ks=[1.0, 1.0])

    # Cannot set the rate constant for a reaction that doesn't exist in the set
    with pytest.raises(ex.CannotSetAttribute):
        rxn_set.set_rate_constant('X', 'Y', k=1.0)
