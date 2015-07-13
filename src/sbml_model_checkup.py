"""Debugging/cleanup checks."""

import itertools
import chen_2009_original_sbml

model = chen_2009_original_sbml.get_model()

# determine redundant species (same name)
sa = sorted(model.species, key=lambda s: (s.label, s.compartment))
rs = [x
      for x in ((n, map(lambda s: s.compartment, it))
          for n, it in itertools.groupby(sa, lambda s: s.label))
      if len(x[1]) > 1 ]

# find reactions where label is not reactants -> products
def format_slist(species):
    return ' + '.join(s.label.replace('endo|', '') for s in species)
def format_rxn(r):
    return ' -> '.join(map(format_slist, [r.reactants, r.products]))
mislabeled_rxns = [(r, fr) for r, fr in
                   zip(model.reactions, map(format_rxn, model.reactions))
                   if r.label != fr and not
                   (' + ->' in r.label and len(r.reactants) == 1)]


print
print "Mismatched reactions:"
print
for r, fr in mislabeled_rxns:
    print "    label:    {}".format(r.label)
    print "    internal: {}".format(fr)
    print

print "Duplicated species:"
for label, compartments in rs:
    print "    {} @ {}".format(label, ', '.join(compartments))
