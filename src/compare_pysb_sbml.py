"""Print SBML model species and PySB model species for 1-1 comparison."""

import difflib
import re
import chen_2009_original_sbml
import chen_2009
from pysb.bng import generate_equations

def get_pysb_names():
    """Return species names based on original SBML model naming convention."""

    model = chen_2009.model
    generate_equations(model)

    ordering = ('EGF HRG ErbB1 ErbB2 ErbB3 ErbB4 ATP RTK_Pase GAP Shc Grb2 Gab1 Shp2 '
                'Pase9t PI3K PIP2 Sos MEK Pase2 ERK Pase3 Raf Pase1 Ras cPP PIP3 '
                'AKT Pase4 PDK1 PTEN Shp').split()
    ordering_map = {p: i for i, p in enumerate(ordering)}

    all_species = [str(i) for i in model.species]
    all_species = [i for i in all_species  if i not in ('__source()', '__sink()')]

    monomers_in_species = [i.split(" % ") for i in all_species]

    species_list = []
    for mlist in monomers_in_species:
         mlist = [s.replace('Erb', 'ErbB') for s in mlist]
         mlist = [s.replace('SHC', 'Shc') for s in mlist]
         mlist = [s.replace('SOS', 'Sos') for s in mlist]
         mlist = [s.replace('Pase_9t', 'Pase9t') for s in mlist]
         mlist = [s.replace('RTK', 'RTK_Pase') for s in mlist]
         mlist = [s.replace('RAF', 'Raf') for s in mlist]
         mlist = [s.replace('RAS', 'Ras') for s in mlist]
         temp = sorted(mlist, key=lambda s: ordering_map[s[:s.index('(')]])
         temp = [re.sub(r'([^(]+).*state=\'p\'.*', r'\1#P', i) for i in temp]
         temp = [re.sub(r'([^(]+).*state=\'pp\'.*', r'\1#P#P', i) for i in temp]
         temp = [re.sub(r'([^(]+).*state=\'gdp\'.*', r'\1:GDP', i) for i in temp]
         temp = [re.sub(r'([^(]+).*state=\'gtp\'.*', r'\1:GTP', i) for i in temp]
         temp = [re.sub(r'([^(]+).*state=\'active_gtp\'.*', r'\1_activated:GTP', i) for i in temp]
         temp = [re.sub(r'(Raf).*state=\'p_ser\'.*', r'\1:P:Ser', i) for i in temp]
         temp = [re.sub(r'([^(]+).*', r'\1', i) for i in temp]
         s = ':'.join(temp)
         s = re.sub(r'ATP:(GAP:Grb2:Gab1)', r'\1:ATP', s)
         s = re.sub(r'EGF:EGF', r'EGF', s)
         s = re.sub(r'(PIP2:?){2,}', lambda m: '(PIP2)%d' % ((len(m.group())+1)/5), s)
         s = re.sub(r'(ErbB\d)#P:(ErbB\d)#P', r'(\1:\2)#P', s)
         s = re.sub(r'^\(ErbB2:ErbB([34])\)', r'(ErbB\1:ErbB2)', s)
         s = re.sub(r'^\(ErbB2:ErbB2\)', r'2(ErbB2)', s)
         s = re.sub(r'EGF:\(ErbB1:ErbB1\)', r'2(EGF:ErbB1)', s)
         if re.match(r'2\(EGF:ErbB1\).*Gab1#[^:]+(?!:PI3K$|:PI3K:PIP2$)', s):
              s = re.sub(r'(Gab1#[^:]+)', r'(\1)', s)
         s = re.sub(r'(Shc#P)', r'(\1)', s)
         s = re.sub(r'(Sos:)(Ras:G[DT]P)', r'\1(\2)', s)
         s = re.sub(r'AKT#P#P', r'AKT:P:P', s)
         s = re.sub(r'HRG:ErbB1:ErbB([34])', r'(HRG:ErbB\1:ErbB1)', s)
         s = re.sub(r'HRG:ErbB2:ErbB([34])', r'(HRG:ErbB\1):ErbB2', s)
         species_list.append(s)

    for i, comp in enumerate(all_species):
        if "comp='endo'" in comp:
            species_list[i] = 'endo|' + species_list[i]

    species_list.sort()

    return species_list


def get_sbml_names_species():
    model = chen_2009_original_sbml.get_model()
    ignore = ('_i', '_h', 'Inh')
    species = sorted([s for s in model.species
                      if not any(i in s.label for i in ignore)],
                     key=lambda x: x.label)
    names = [s.label for s in species]
    return names, species


pysb_names = get_pysb_names()
sbml_names, sbml_species = get_sbml_names_species()

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
