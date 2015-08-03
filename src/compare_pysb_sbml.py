"""Print SBML model species and PySB model species for 1-1 comparison."""

from __future__ import division
import difflib
import re
import chen_2009_original_sbml
import chen_2009
from pysb.bng import generate_equations

def get_pysb_species():
    """Return species names mostly aligned to SBML model naming convention."""

    # Note that when we refer to to "sbml model" or just "sbml" below we
    # specifically mean the sbml version of the Chen 2009 model (as opposed to
    # the pysb version of it) and not any sbml model in general.

    model = chen_2009.model
    generate_equations(model)

    # This is a rough total ordering of the protein names as used in the sbml
    # model species labels. E.g. EGF always comes before ErbB1. The sbml naming
    # is not totally consistent which is why this is a "rough" ordering. We will
    # use this ordering below to sort pysb monomer names within each species to
    # produce names ordered as similarly as possible to the sbml ones.
    ordering = ('EGF HRG ErbB1 ErbB2 ErbB3 ErbB4 ATP RTK_Pase GAP Shc Grb2 '
                'Gab1 Shp2 Pase9t PI3K PIP2 Sos ERK MEK Pase2 Pase3 Raf Pase1 '
                'Ras cPP PIP3 AKT Pase4 PDK1 PTEN Shp').split()
    ordering_map = {p: i for i, p in enumerate(ordering)}

    # Augment species objects with their original numbering.
    for i, s in enumerate(model.species):
        s.index = i
    # Throw out source and sink species and convert to strings.
    species = [s for s in model.species
               if str(s) not in ('__source()', '__sink()')]
    species_str = [str(s) for s in species]
    # Split complexes on % to produce monomers.
    monomers_in_species = [i.split(" % ") for i in species_str]

    names = []
    for mlist in monomers_in_species:
        # Fix some minor spelling/case differences in some proteins.
        mlist = [s.replace('Erb', 'ErbB') for s in mlist]
        mlist = [s.replace('SHC', 'Shc') for s in mlist]
        mlist = [s.replace('SOS', 'Sos') for s in mlist]
        mlist = [s.replace('Pase_9t', 'Pase9t') for s in mlist]
        mlist = [s.replace('RTK', 'RTK_Pase') for s in mlist]
        mlist = [s.replace('RAF', 'Raf') for s in mlist]
        mlist = [s.replace('RAS', 'Ras') for s in mlist]
        # Topological sort on monomers based on protein ordering defined above.
        temp = sorted(mlist, key=lambda s: ordering_map[s[:s.index('(')]])
        # Convert various state flags to text suffixes used in sbml.
        temp = [re.sub(r'([^(]+).*state=\'p\'.*', r'\1#P', i) for i in temp]
        temp = [re.sub(r'([^(]+).*state=\'pp\'.*', r'\1#P#P', i) for i in temp]
        temp = [re.sub(r'([^(]+).*state=\'gdp\'.*', r'\1:GDP', i) for i in temp]
        temp = [re.sub(r'([^(]+).*state=\'gtp\'.*', r'\1:GTP', i) for i in temp]
        temp = [re.sub(r'([^(]+).*state=\'active_gtp\'.*', r'\1_activated:GTP', i) for i in temp]
        temp = [re.sub(r'(Raf).*state=\'p_ser\'.*', r'\1:P:Ser', i) for i in temp]
        # Strip remaining sites and parens, leaving just name and suffix.
        temp = [re.sub(r'([^(]+).*', r'\1', i) for i in temp]
        # Join the proteins back together on : as used in sbml.
        s = ':'.join(temp)
        # Apply some special case naming patterns to match sbml.
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
            # For species in the endo compartment, prepend endo prefix to name
            # and apply one special case name fixup.
            names[i] = re.sub(r'(EGF:ErbB1:ErbB[234])$', r'(\1)', names[i])
            names[i] = 'endo|' + names[i]

    # Sort names and species by names.
    names, species = zip(*sorted(zip(names, species)))

    return names, species, model


def get_sbml_species():
    model = chen_2009_original_sbml.get_model()
    # We will ignore species whose labels contain these strings.
    ignore_patterns = ('_i', '_h', 'Inh')
    # We will ignore these individually named species.
    ignore_names = (
        # Degradation sinks.
        'c13', 'c520', 'c86',
        # MEK#P#P:ERK "_i" species that are missing _i in label.
        'c80', 'c82',
        )
    species = [s for s in model.species
               if not any(i in s.label for i in ignore_patterns)
               and s.name not in ignore_names]
    names = [s.label for s in species]

    # Sort names and species by names.
    names, species = zip(*sorted(zip(names, species)))

    return names, species, model


pysb_names, pysb_species, pysb_model = get_pysb_species()
sbml_names, sbml_species, sbml_model = get_sbml_species()

if len(pysb_names) != len(set(pysb_names)):
    raise RuntimeError("Duplicate pysb species names")
if len(sbml_names) != len(set(sbml_names)):
    raise RuntimeError("Duplicate sbml species names")

# Whether or not to print matching lines (mismatches will always print).
print_matches = False

# Count matches.
species_matches = reaction_matches = 0

fmt = '%-60s\t%1s\t%-51s'
print fmt % ('SBML', '', 'PySB')
print fmt % ('=' * 50, '', '=' * 50)
print
# Use difflib to calculate the difference between our name lists.
sm = difflib.SequenceMatcher(None, sbml_names, pysb_names)
for tag, i1, i2, j1, j2 in sm.get_opcodes():
    if tag in ('delete', 'replace'):
        # Species in sbml model but not pysb model.
        for si in range(i1, i2):
            ss = '%s - %s' % (sbml_names[si], sbml_species[si].name)
            print fmt % (ss, '', '')
    if tag in ('insert', 'replace'):
        # Species in pysb model but not in sbml model.
        for sj in range(j1, j2):
            ss = '%s - s%d' % (pysb_names[sj], sj)
            print fmt % ('', '', ss)
    if tag == 'equal':
        # Species in both.
        species_matches += i2 - i1
        if print_matches:
            for si, sj in zip(range(i1, i2), range(j1, j2)):
                sbml = '%s - %s' % (sbml_names[si], sbml_species[si].name)
                pysb = '%s - s%d' % (pysb_names[sj], sj)
                print fmt % (sbml, '=', pysb)

pysb_to_sbml = {p.index: s.name for p, s in zip(pysb_species, sbml_species)}
pysb_to_sbml[46] = '(degraded)'
pysb_reactions = sorted(
    ' -> '.join(sorted([
            ' + '.join(sorted(pysb_to_sbml[s] for s in r[field]))
            for field in 'reactants', 'products'
            ], key=lambda x: 999 if x == '(degraded)' else -len(x)))
    for r in chen_2009.model.reactions_bidirectional)
sinks = [s for s in sbml_model.species if s.name in ('c13', 'c520', 'c86')]
sbml_reactions = sorted(
    ' -> '.join([
            ' + '.join(sorted(s.name if s not in sinks else '(degraded)' for s in getattr(r, field)))
            for field in 'reactants', 'products'
            ])
    for r in sbml_model.reactions
    if all(s in sbml_species or s in sinks for s in r.reactants + r.products))
sm = difflib.SequenceMatcher(None, sbml_reactions, pysb_reactions)
for tag, i1, i2, j1, j2 in sm.get_opcodes():
    if tag in ('delete', 'replace'):
        # Reactions in sbml model but not pysb model.
        for si in range(i1, i2):
            ss = '%s - %s' % (sbml_reactions[si], sbml_model.reactions[si].name)
            print fmt % (ss, '', '')
    if tag in ('insert', 'replace'):
        # Reactions in pysb model but not in sbml model.
        for sj in range(j1, j2):
            ss = '%s - r%d' % (pysb_reactions[sj], sj)
            print fmt % ('', '', ss)
    if tag == 'equal':
        # Reactions in both.
        reaction_matches += i2 - i1

species_match_percent = species_matches / len(sbml_names) * 100
print
print "Species matches: %d / %d -- %.2f%% %s" % (
    species_matches, len(sbml_names), species_match_percent,
    u'\U0001f37b' if species_match_percent == 100 else ''
    )
print "SBML species missed: %d" % (len(sbml_names) - species_matches)
print "PySB surplus species: %d" % (len(pysb_names) - species_matches)

reaction_match_percent = reaction_matches / len(sbml_reactions) * 100
print
print "Reaction matches: %d / %d -- %.2f%% %s" % (
    reaction_matches, len(sbml_reactions), reaction_match_percent,
    u'\U0001f37b' if reaction_match_percent == 100 else ''
    )
print "SBML reactions missed: %d" % (len(sbml_reactions) - reaction_matches)
print "PySB surplus reactions: %d" % (len(pysb_reactions) - reaction_matches)
