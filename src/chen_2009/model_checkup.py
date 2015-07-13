"""Debugging/cleanup checks."""

import itertools
import original

model = original.get_model()

# determine truly redundant species -- same name, same compartment
sa = sorted(model.species, key=lambda s: (s.label, s.compartment))
rs = [x
      for x in ((n, map(lambda s: s.compartment, it))
          for n, it in itertools.groupby(sa, lambda s: s.label))
      if len(x[1]) == 2 and x[1][0] == x[1][1]]

# find reactions where product is not a trivial concatenation of reactants
# (e.g. A + B -> A:B)
mismatch_rxns = [r for r in model.reactions
                 if r.products[0].label != ':'.join([s.label for s in r.reactants])]
