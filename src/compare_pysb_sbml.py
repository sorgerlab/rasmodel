"""Print SBML model species and PySB model species for 1-1 comparison."""

import difflib
import chen_2009_original_sbml
import chen_2009

ignore = ('_i', '_h', 'Inh')
sbml_species = sorted([s.label for s in
                       chen_2009_original_sbml.get_model().species
                       if not any(i in s.label for i in ignore)])
pysb_species = sorted(chen_2009.get_species_names())

fmt = '%-51s\t%-51s'
print fmt % ('SBML', 'PySB')
print fmt % (('=' * 51,) * 2)
print
sm = difflib.SequenceMatcher(None, sbml_species, pysb_species)
for tag, i1, i2, j1, j2 in sm.get_opcodes():
    if tag == 'delete':
        for s in sbml_species[i1:i2]:
            print fmt % (s, '')
    elif tag == 'insert':
        for s in pysb_species[j1:j2]:
            print fmt % ('', s)
