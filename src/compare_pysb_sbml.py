"""Print SBML model species and PySB model species for 1-1 comparison."""

from __future__ import division
import difflib
import re
import argparse
import sys
import chen_2009_original_sbml
from REM import chen_2009
import pysb
from pysb.bng import generate_equations

def get_pysb_species():
    """Return species names mostly aligned to SBML model naming convention."""

    # Note that when we refer to to "sbml model" or just "sbml" below we
    # specifically mean the sbml version of the Chen 2009 model (as opposed to
    # the pysb version of it) and not any sbml model in general.

    model = chen_2009.model

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

    labels = []
    for mlist in monomers_in_species:
        # Fix some minor spelling/case differences in some proteins.
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
        temp = [re.sub(r'([^(]+).*state=\'full_act\'.*', r'\1-FullActive', i) for i in temp]
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
        if '-FullActive' in s:
            s = s.replace('-FullActive', '')
            s = s + '-FullActive'
        s = re.sub(r'^EGF:ErbB1:ErbB1:ATP:ATP(-FullActive|)', r'2(EGF:ErbB1:ATP)\1', s)
        labels.append(s)

    for i, comp in enumerate(species_str):
        if "comp='endo'" in comp:
            # For species in the endo compartment, prepend endo prefix to name
            # and apply one special case name fixup.
            labels[i] = re.sub(r'(EGF:ErbB1:ErbB[234])$', r'(\1)', labels[i])
            labels[i] = 'endo|' + labels[i]

    ics = [''] * len(species)
    for ic_species, ic_parameter in model.initial_conditions:
        if ic_parameter.value != 0 and str(ic_species) != '__source()':
            idx = next(i for i, s in enumerate(species)
                       if s.is_equivalent_to(ic_species))
            ics[idx] = ' @ %.17g' % ic_parameter.value

    names = [label + ic for ic, label in zip(ics, labels)]

    # Sort names and species by names.
    names, species = zip(*sorted(zip(names, species)))

    return names, species


def get_sbml_species():
    model = chen_2009_original_sbml.model
    # We will ignore species whose labels contain these strings.
    ignore_patterns = ('_i', '_h', 'Inh')
    # We will ignore these individually named species.
    ignore_names = (
        # Degradation sinks.
        'c13', 'c520', 'c86',
        # MEK#P#P:ERK "_i" species that are missing _i in label.
        'c80', 'c82', 'c96', 'c98',
        )
    species = [s for s in model.species
               if not any(i in s.label for i in ignore_patterns)
               and s.name not in ignore_names]
    labels = [s.label for s in species]
    ics = [' @ %.17g' % s.initial_amount if s.initial_amount != 0 else ''
           for s in species]
    names = [label + ic for label, ic in zip(labels, ics)]

    # Sort names and species by names.
    names, species = zip(*sorted(zip(names, species)))

    return names, species


def get_pysb_reactions():
    model = chen_2009.model

    pysb_to_sbml = {p.index: s.name for p, s in zip(pysb_species, sbml_species)}
    sink_index = next(i for i, s in enumerate(model.species)
                      if str(s) == '__sink()')
    pysb_to_sbml[sink_index] = '(degraded)'

    reactions = model.reactions_bidirectional
    for i, r in enumerate(reactions):
        r['index'] = i

    def format_side(indexes):
        labels = [pysb_to_sbml[i] for i in indexes]
        return ' + '.join(sorted(labels))

    def format_param(parameter):
        if parameter is not None and parameter.get_value() != 0:
            pname = parameter.name
            # Fixup for BNG symmetry corrections.
            if pname.endswith('_symmetric'):
                orig_param = parameter.expr / 2
                if not isinstance(orig_param, pysb.Parameter):
                    raise RuntimeError("Unexpected expression structure")
                pname = orig_param.name
            return pname
        else:
            return '<0>'

    descriptions = []
    for reaction in reactions:
        assert len(reaction['rule']) == 1
        rule = model.rules[reaction['rule'][0]]
        left = reaction['reactants']
        right = reaction['products']
        kf = rule.rate_forward
        kr = rule.rate_reverse
        # Swap sides to put smaller species list or '(degraded)' on the right.
        if len(left) < len(right) or left == (sink_index,):
            left, right = right, left
            kf, kr = kr, kf
        desc = '%s -> %s {%s, %s}' % (format_side(left), format_side(right),
                                      format_param(kf), format_param(kr))
        descriptions.append(desc)

    # Sort names and reactions by names.
    descriptions, reactions = zip(*sorted(zip(descriptions, reactions)))

    return descriptions, reactions


def get_sbml_reactions():
    model = chen_2009_original_sbml.model

    sinks = tuple(s for s in model.species if s.name in ('c13', 'c520', 'c86'))

    _, sbml_species = get_sbml_species()
    wanted = sbml_species + sinks
    # Skip reactions that aren't fully in our scope.
    reactions = [r for r in model.reactions
                 if all(s in wanted for s in r.reactants + r.products)]

    def format_side(species):
        labels = [s.name if s not in sinks else '(degraded)' for s in species]
        return ' + '.join(sorted(labels))

    def format_param(parameter):
        return parameter.name if parameter.value != 0 else '<0>'

    descriptions = []
    for reaction in reactions:
        left = reaction.reactants
        right = reaction.products
        kf = reaction.kf
        kr = reaction.kr
        desc = '%s -> %s {%s, %s}' % (format_side(left), format_side(right),
                                      format_param(kf), format_param(kr))
        descriptions.append(desc)

    # Sort names and reactions by names.
    descriptions, reactions = zip(*sorted(zip(descriptions, reactions)))

    return descriptions, reactions


argparser = argparse.ArgumentParser()
argparser.add_argument('-m', '--print-matches', action='store_true',
                       help="Print matching elements too (mismatches are "
                       "always printed)")
args = argparser.parse_args(sys.argv[1:])

pysb_model = chen_2009.model
generate_equations(pysb_model)

chen_2009_original_sbml.load_model()
sbml_model = chen_2009_original_sbml.model

pysb_species_names, pysb_species = get_pysb_species()
sbml_species_names, sbml_species = get_sbml_species()

if len(pysb_species_names) != len(set(pysb_species_names)):
    raise RuntimeError("Duplicate pysb species names")
if len(sbml_species_names) != len(set(sbml_species_names)):
    raise RuntimeError("Duplicate sbml species names")

# Count matches.
species_matches = reaction_matches = 0

fmt = '%-60s\t%1s\t%-51s'
print fmt % ('SBML', '', 'PySB')
print fmt % ('=' * 50, '', '=' * 50)
print
# Use difflib to calculate the difference between our name lists.
sm = difflib.SequenceMatcher(None, sbml_species_names, pysb_species_names)
for tag, i1, i2, j1, j2 in sm.get_opcodes():
    if tag in ('delete', 'replace'):
        # Species in sbml model but not pysb model.
        for si in range(i1, i2):
            ss = '%s : %s' % (sbml_species_names[si], sbml_species[si].name)
            print fmt % (ss, '', '')
    if tag in ('insert', 'replace'):
        # Species in pysb model but not in sbml model.
        for sj in range(j1, j2):
            ss = '%s : s%d' % (pysb_species_names[sj], sj)
            print fmt % ('', '', ss)
    if tag == 'equal':
        # Species in both.
        species_matches += i2 - i1
        if args.print_matches:
            for si, sj in zip(range(i1, i2), range(j1, j2)):
                sbml = '%s : %s' % (sbml_species_names[si], sbml_species[si].name)
                pysb = '%s : s%d' % (pysb_species_names[sj], sj)
                print fmt % (sbml, '=', pysb)
species_match_percent = species_matches / len(sbml_species_names) * 100

if species_match_percent == 100.0:
    pysb_reaction_descs, pysb_reactions = get_pysb_reactions()
    sbml_reaction_descs, sbml_reactions = get_sbml_reactions()
    sm = difflib.SequenceMatcher(None, sbml_reaction_descs, pysb_reaction_descs)
    for tag, i1, i2, j1, j2 in sm.get_opcodes():
        if tag in ('delete', 'replace'):
            # Reactions in sbml model but not pysb model.
            for si in range(i1, i2):
                ss = '%s : %s' % (sbml_reaction_descs[si],
                                  sbml_reactions[si].name)
                print fmt % (ss, '', '')
        if tag in ('insert', 'replace'):
            # Reactions in pysb model but not in sbml model.
            for sj in range(j1, j2):
                ss = '%s : r%d, %s' % (pysb_reaction_descs[sj],
                                       pysb_reactions[sj]['index'],
                                       pysb_reactions[sj]['rule'][0])
                print fmt % ('', '', ss)
        if tag == 'equal':
            # Reactions in both.
            reaction_matches += i2 - i1
            if args.print_matches:
                for si, sj in zip(range(i1, i2), range(j1, j2)):
                    sbml = '%s : %s' % (sbml_reaction_descs[si],
                                        sbml_reactions[si].name)
                    pysb = '%s : r%d, %s' % (pysb_reaction_descs[sj],
                                             pysb_reactions[sj]['index'],
                                             pysb_reactions[sj]['rule'][0])
                    print fmt % (sbml, '=', pysb)
    reaction_match_percent = reaction_matches / len(sbml_reactions) * 100
else:
    reaction_match_percent = 0

if reaction_match_percent == 100.0:
    pysb_parameters = pysb_model.parameters_rules()
    sbml_parameters = {}
    sbml_model.export_globals(sbml_parameters)
    parameter_mismatches = []
    for pysb_parameter in pysb_parameters:
        # Fixup for BNG symmetry corrections.
        pname = pysb_parameter.name.replace('_symmetric', '')
        try:
            sbml_parameter = sbml_parameters[pname]
        except KeyError:
            parameter_mismatches.append((pysb_parameter, None))
            continue
        if pysb_parameter.get_value() != sbml_parameter.value:
            parameter_mismatches.append((pysb_parameter, sbml_parameter))
    num_param_matches = len(pysb_parameters) - len(parameter_mismatches)
    parameter_match_percent = num_param_matches / len(pysb_parameters) * 100

print
print "Species matches: %d / %d -- %.2f%% %s" % (
    species_matches, len(sbml_species), species_match_percent,
    u'\U0001f37b' if species_match_percent == 100 else ''
    )
print "SBML species missed: %d" % (len(sbml_species) - species_matches)
print "PySB surplus species: %d" % (len(pysb_species) - species_matches)

if species_match_percent == 100.0:
    print
    print "Reaction matches: %d / %d -- %.2f%% %s" % (
        reaction_matches, len(sbml_reactions), reaction_match_percent,
        u'\U0001f37b' if reaction_match_percent == 100 else ''
        )
    print "SBML reactions missed: %d" % (len(sbml_reactions) - reaction_matches)
    print "PySB surplus reactions: %d" % (len(pysb_reactions) - reaction_matches)
else:
    print "\nSkipping reaction comparison until species match is 100%"

if reaction_match_percent == 100.0:
    print
    print "Parameter value matches: %d / %d -- %.2f%% %s" % (
        num_param_matches, len(pysb_parameters), parameter_match_percent,
        u'\U0001f37b' if parameter_match_percent == 100 else ''
        )

# I should probably have generated this list by running the regex from the
# original simulation script against the original species names, but deadlines
# are approaching! -JLM
# Note that this list has already been converted to PySB model species indexes!
sbml_perbb1_species = set([
    85, 86, 110, 111, 119, 120, 121, 122, 153, 154, 156, 157, 158, 159, 160,
    161, 162, 163, 166, 167, 169, 170, 172, 173, 175, 176, 188, 189, 202, 205,
    206, 207, 208, 209, 210, 211, 212, 213, 216, 217, 218, 219, 220, 221, 222,
    223, 224, 225, 226, 227, 228, 229, 231, 232, 234, 235, 236, 237, 238, 239,
    240, 241, 246, 247, 251, 252, 253, 254, 255, 256, 257, 258, 263, 264, 270,
    271, 272, 273, 274, 275, 279, 280, 281, 282, 283, 284, 285, 286, 287, 288,
    289, 290, 291, 292, 293, 294, 295, 296, 297, 298, 299, 300, 307, 308, 309,
    310, 311, 312, 313, 314, 315, 316, 317, 318, 319, 320, 321, 322, 327, 328,
    329, 330, 337, 338, 339, 340, 344, 345, 346, 347, 348, 349, 350, 351, 352,
    353, 354, 355, 356, 357, 358, 359, 368, 369, 370, 371, 373, 374, 379, 380,
    381, 382, 383, 384, 387, 388, 392, 393, 395, 396, 408, 409, 411, 412, 414,
    415, 419, 420, 425, 426, 427, 428, 431, 432
    ])
# Should also be comparing coefficients to original model but again, no time...
pysb_perbb1_species = (
    set(pysb_model.observables['pErbB1_total'].species) -
    set(pysb_model.observables['pErbB11_exceptions'].species) -
    set(pysb_model.observables['pErbB12_exceptions'].species)
    )
print
if pysb_perbb1_species != sbml_perbb1_species:
    print "pErbB1 observable is WRONG:"
    print "  missing:", ', '.join(str(i) for i in sbml_perbb1_species - pysb_perbb1_species)
    print "  extra:", ', '.join(str(i) for i in pysb_perbb1_species - sbml_perbb1_species)
else:
    print "pErbB1 observable is correct"

# For debugging/testing.
sbml2pysb = {x[0].name: x[1].index
             for x in sorted(zip(sbml_species, pysb_species),
                             key=lambda x: int(x[0].name[1:]))}
pysb2sbml = {x[1].index: x[0].name
             for x in sorted(zip(sbml_species, pysb_species),
                             key=lambda x: x[1].index)}
