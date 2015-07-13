"""Print SBML model species and PySB model species for 1-1 comparison."""

import difflib
import chen_2009_original_sbml
import chen_2009

ignore = ('_i', '_h', 'Inh')
sbml_species = sorted([s for s in chen_2009_original_sbml.get_model().species
                       if not any(i in s.label for i in ignore)],
                      key=lambda x: x.label)
sbml_names = [s.label for s in sbml_species]
pysb_names = sorted(chen_2009.get_species_names())

fmt = '%-60s\t%-51s'
print fmt % ('SBML', 'PySB')
print fmt % (('=' * 50,) * 2)
print
sm = difflib.SequenceMatcher(None, sbml_names, pysb_names)
for tag, i1, i2, j1, j2 in sm.get_opcodes():
    if tag == 'delete':
        for si in xrange(i1, i2):
            ss = '%s - %s' % (sbml_names[si], sbml_species[si].name)
            print fmt % (ss, '')
    elif tag == 'insert':
        for s in pysb_names[j1:j2]:
            print fmt % ('', s)
