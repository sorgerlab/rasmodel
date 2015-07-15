"""Print SBML model species and PySB model species for 1-1 comparison."""

from __future__ import division
import difflib
import re
import chen_2009_original_sbml
import chen_2009
from pysb.bng import generate_equations

def get_pysb_species():
    """Return species names mostly aligned to SBML model naming convention."""

    model = chen_2009.model
    generate_equations(model)

    ordering = ('EGF HRG ErbB1 ErbB2 ErbB3 ErbB4 ATP RTK_Pase GAP Shc Grb2 '
                'Gab1 Shp2 Pase9t PI3K PIP2 Sos ERK MEK Pase2 Pase3 Raf Pase1 '
                'Ras cPP PIP3 AKT Pase4 PDK1 PTEN Shp').split()
    ordering_map = {p: i for i, p in enumerate(ordering)}

    species = [s for s in model.species
               if str(s) not in ('__source()', '__sink()')]
    species_str = [str(s) for s in species]

    monomers_in_species = [i.split(" % ") for i in species_str]

    names = []
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
         s = re.sub(r'^EGF:ErbB1:ErbB1:ATP$', r'2(EGF:ErbB1:ATP)', s)
         s = re.sub(r'(Shc#P)', r'(\1)', s)
         s = re.sub(r'(Sos:)(Ras:G[DT]P)', r'\1(\2)', s)
         s = re.sub(r'AKT#P#P', r'AKT:P:P', s)
         s = re.sub(r'HRG:ErbB1:ErbB([34])', r'(HRG:ErbB\1:ErbB1)', s)
         s = re.sub(r'HRG:ErbB2:ErbB([34])', r'(HRG:ErbB\1):ErbB2', s)
         s = re.sub(r'^EGF:ErbB1:ErbB1:ATP:ATP$', r'2(EGF:ErbB1:ATP)-FullActive', s)
         names.append(s)

    for i, comp in enumerate(species_str):
        if "comp='endo'" in comp:
            names[i] = re.sub(r'(EGF:ErbB1:ErbB[234])$', r'(\1)', names[i])
            names[i] = 'endo|' + names[i]

    names, species = zip(*sorted(zip(names, species)))

    return names, species


def get_sbml_species():
    model = chen_2009_original_sbml.get_model()
    ignore_patterns = ('_i', '_h', 'Inh')
    ignore_names = (
        'c13', 'c520', 'c86', # Degradation sinks
        'c80', 'c82', # MEK#P#P:ERK "_i" species that are missing _i in label
        )
    species = [s for s in model.species
               if not any(i in s.label for i in ignore_patterns)
               and s.name not in ignore_names]
    names = [s.label for s in species]
    names, species = zip(*sorted(zip(names, species)))
    return names, species


pysb_names, pysb_species = get_pysb_species()
sbml_names, sbml_species = get_sbml_species()

if len(pysb_names) != len(set(pysb_names)):
    raise RuntimeError("Duplicate pysb species names")
if len(sbml_names) != len(set(sbml_names)):
    raise RuntimeError("Duplicate sbml species names")

print_matches = False

matches = 0
fmt = '%-60s\t%1s\t%-51s'
print fmt % ('SBML', '', 'PySB')
print fmt % ('=' * 50, '', '=' * 50)
print
sm = difflib.SequenceMatcher(None, sbml_names, pysb_names)
for tag, i1, i2, j1, j2 in sm.get_opcodes():
    if tag in ('delete', 'replace'):
        for si in range(i1, i2):
            ss = '%s - %s' % (sbml_names[si], sbml_species[si].name)
            print fmt % (ss, '', '')
    if tag in ('insert', 'replace'):
        for sj in range(j1, j2):
            ss = '%s - s%d' % (pysb_names[sj], sj)
            print fmt % ('', '', ss)
    if tag == 'equal':
        matches += i2 - i1
        if print_matches:
            for si, sj in zip(range(i1, i2), range(j1, j2)):
                sbml = '%s - %s' % (sbml_names[si], sbml_species[si].name)
                pysb = '%s - s%d' % (pysb_names[sj], sj)
                print fmt % (sbml, '=', pysb)

print
print "Total matches: %d / %d (%.2f%%)" % (matches, len(sbml_names),
                                         matches / len(sbml_names) * 100)
print "SBML species missed: %d" % (len(sbml_names) - matches)
print "PySB surplus species: %d" % (len(pysb_names) - matches)
