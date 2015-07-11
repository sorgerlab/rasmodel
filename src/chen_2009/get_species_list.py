from pysb.bng import generate_equations
from chen_2009 import model
import re

generate_equations(model)

ordering = ('EGF HRG ErbB1 ErbB2 ErbB3 ErbB4 ATP RTK GAP SHC Grb2 Gab1 Shp2 '
            'Pase_9t PI3K PIP2 SOS MEK Pase2 ERK Pase3 RAF Pase1 RAS cPP PIP3 '
            'AKT Pase4 PDK1 PTEN Shp').split()
ordering_map = {p: i for i, p in enumerate(ordering)}

all_species = [str(i) for i in model.species]

monomers_in_species = [i.split(" % ") for i in all_species]

species_list = []
for mlist in monomers_in_species:
     if mlist in (['__source()'], ['__sink()']):
          continue
     mlist = [s.replace('Erb', 'ErbB') for s in mlist]
     temp = sorted(mlist, key=lambda s: ordering_map[s[:s.index('(')]])
     temp = [re.sub(r'([^(]+).*state=\'p\'.*', r'\1#P', i) for i in temp]
     temp = [re.sub(r'([^(]+).*state=\'pp\'.*', r'\1#P#P', i) for i in temp]
     temp = [re.sub(r'([^(]+).*state=\'gdp\'.*', r'\1:GDP', i) for i in temp]
     temp = [re.sub(r'([^(]+).*state=\'gtp\'.*', r'\1:GTP', i) for i in temp]
     temp = [re.sub(r'([^(]+).*state=\'active_gtp\'.*', r'(\1_activated:GTP)', i) for i in temp]
     temp = [re.sub(r'(RAF).*state=\'p_ser\'.*', r'\1:P:Ser', i) for i in temp]
     temp = [re.sub(r'([^(]+).*', r'\1', i) for i in temp]
     s = ':'.join(temp)
     s = re.sub(r'EGF:EGF', r'EGF', s)
     s = re.sub(r'(PIP2:?){2,}', lambda m: '(PIP2)%d' % ((len(m.group())+1)/5), s)
     s = re.sub(r'(ErbB\d)#P:(ErbB\d)#P', r'(\1:\2)#P', s)
     species_list.append(s)

for i, comp in enumerate(all_species):
    if "comp='endo'" in comp:
        species_list[i] = 'endo|' + species_list[i]

print '\n'.join(species_list)
